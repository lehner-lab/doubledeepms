#!/usr/bin/env python

#######################################################################
## COMMANDLINE ARGUMENTS ##
#######################################################################

import argparse

#Create parser
parser = argparse.ArgumentParser()
 
#Add arguments to the parser
parser.add_argument("--data_train", help = "Training data")
parser.add_argument("--data_valid", default = None, help = "Validation data")
parser.add_argument("--data_obs", help = "All observed data")
parser.add_argument("--output_directory", "-o")
parser.add_argument("--number_additive_traits", "-n", default = 1, type = int, help = "Number of additive traits")
parser.add_argument("--l1_regularization_factor", default = "0.0001,0.001,0.01,0.1", help = "L1 regularization factor for binding additive trait layer (default:0.0001,0.001,0.01,0.1)")
parser.add_argument("--l2_regularization_factor", default = "0.0001,0.001,0.01,0.1", help = "L2 regularization factor for binding additive trait layer (default:0.0001,0.001,0.01,0.1)")
parser.add_argument("--num_epochs_grid", "-e", default = 100, type = int, help = "Number of epochs to train the model during grid search")
parser.add_argument("--num_epochs", "-p", default = 10000, type = int, help = "Maximum number of epochs to train the final model")
parser.add_argument("--num_samples", "-s", default = "128,256,512,1024", help = "Number of samples per gradient update (default:128,256,512,1024)")
parser.add_argument("--learning_rate", "-a", default = "0.0001,0.001,0.01,0.1", help = "Learning rate (default:0.0001,0.001,0.01,0.1)")
parser.add_argument("--num_resamplings", "-r", default = 10, type = int, help = "Number of random resamples from fitness distribution (default:10)")
parser.add_argument("--early_stopping", "-l", default = False, type = bool, help = "Whether to stop early (default:False)")
parser.add_argument("--num_models", "-u", default = 10, type = int, help = "Number of final models to fit (default:10)")
parser.add_argument("--random_seed", "-d", default = 1, type = int, help = "Random seed (default:1)")

#Parse the arguments
args = parser.parse_args()
print(args)

data_train_file = args.data_train 
data_valid_file = args.data_valid 
data_obs_file = args.data_obs 
output_directory = args.output_directory 
number_additive_traits = args.number_additive_traits
num_epochs_grid = args.num_epochs_grid
num_epochs = args.num_epochs
num_resamplings = args.num_resamplings
early_stopping = args.early_stopping
num_models = args.num_models
random_seed = args.random_seed
#Grid search arguments
l1 = [float(i) for i in args.l1_regularization_factor.split(",")]
l2 = [float(i) for i in args.l2_regularization_factor.split(",")]
batch_size = [int(i) for i in args.num_samples.split(",")]
learn_rate = [float(i) for i in args.learning_rate.split(",")]

#######################################################################
## PACKAGES ##
#######################################################################

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model
import random
from sklearn.metrics import mean_absolute_error
import matplotlib
import os

from tensorflow.keras.constraints import Constraint
from tensorflow.keras.layers import Layer
from tensorflow.keras import backend as K

#######################################################################
## CLASSES ##
#######################################################################

class State_prob_folded(Layer):
    def __init__(self, trainable=False, **kwargs):
        super(State_prob_folded, self).__init__(**kwargs)
        self.supports_masking = True
        self.trainable = trainable
    def build(self, input_shape):
        super(State_prob_folded, self).build(input_shape)
    def call(self, inputs, mask=None):
        return K.pow(K.constant(1.) + K.exp(inputs),-1)
    def get_config(self):
        config = {'trainable': self.trainable}
        base_config = super(State_prob_folded, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))
    def compute_output_shape(self, input_shape):
        return input_shape

class State_prob_bound(Layer):
    def __init__(self, trainable=False, **kwargs):
        super(State_prob_bound, self).__init__(**kwargs)
        self.supports_masking = True
        self.trainable = trainable
    def build(self, input_shape):
        super(State_prob_bound, self).build(input_shape)
    def call(self, inputs, mask=None):
        return K.pow(K.constant(1.) + K.exp(inputs[0]) * (K.constant(1.)+K.exp(inputs[1])),-1)
    def get_config(self):
        config = {'trainable': self.trainable}
        base_config = super(State_prob_bound, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))
    def compute_output_shape(self, input_shape):
        return input_shape

class Between(Constraint):
    def __init__(self, min_value, max_value):
        self.min_value = min_value
        self.max_value = max_value
    def __call__(self, w):        
        return K.clip(w, self.min_value, self.max_value)
    def get_config(self):
        return {'min_value': self.min_value,
                'max_value': self.max_value}

#######################################################################
## FUNCTIONS ##
#######################################################################

#Load model data into sparse tensors
def load_model_data(file_dict):
  data_dict = {}
  for name in file_dict.keys():
    #Initialise
    data_dict[name] = {}
    #Column names
    ALL_COLUMNS = list(pd.read_csv(file_dict[name], nrows = 1).columns)
    SELECT_COLUMNS = [i for i in range(len(ALL_COLUMNS)) if str.startswith(ALL_COLUMNS[i], "dataset_")]
    FOLD_COLUMNS = [i for i in range(len(ALL_COLUMNS)) if str.startswith(ALL_COLUMNS[i], "fold_") or ALL_COLUMNS[i]=="WT"]
    BIND_COLUMNS = [i for i in range(len(ALL_COLUMNS)) if str.startswith(ALL_COLUMNS[i], "bind_") or ALL_COLUMNS[i]=="WT"]
    TARGET_COLUMN = [i for i in range(len(ALL_COLUMNS)) if ALL_COLUMNS[i]=="fitness"]
    TARGET_SD_COLUMN = [i for i in range(len(ALL_COLUMNS)) if ALL_COLUMNS[i]=="fitness_sd"]
    SEQUENCE_COLUMN = [i for i in range(len(ALL_COLUMNS)) if ALL_COLUMNS[i]=="variant_sequence"]
    TRAINING_SET_COLUMN = [i for i in range(len(ALL_COLUMNS)) if ALL_COLUMNS[i]=="training_set"]
    #Save (sparse) tensors
    data_dict[name]["select"] = tf.convert_to_tensor(np.asarray(pd.read_csv(file_dict[name], usecols = SELECT_COLUMNS)), np.float32)
    data_dict[name]["fold"] = tf.sparse.from_dense(tf.convert_to_tensor(np.asarray(pd.read_csv(file_dict[name], usecols = FOLD_COLUMNS)), np.float32))
    data_dict[name]["bind"] = tf.sparse.from_dense(tf.convert_to_tensor(np.asarray(pd.read_csv(file_dict[name], usecols = BIND_COLUMNS)), np.float32))
    data_dict[name]["target"] = tf.convert_to_tensor(np.asarray(pd.read_csv(file_dict[name], usecols = TARGET_COLUMN)), np.float32)
    data_dict[name]["target_sd"] = tf.convert_to_tensor(np.asarray(pd.read_csv(file_dict[name], usecols = TARGET_SD_COLUMN)), np.float32)
    #Save remaining columns
    if len(SEQUENCE_COLUMN)!=0 and len(TRAINING_SET_COLUMN)!=0:
      data_dict[name]["sequence"] = np.asarray(pd.read_csv(file_dict[name], usecols = SEQUENCE_COLUMN))
    if len(TRAINING_SET_COLUMN)!=0:
      data_dict[name]["training_set"] = np.asarray(pd.read_csv(file_dict[name], usecols = TRAINING_SET_COLUMN))
    data_dict[name]["fold_colnames"] = np.asarray([ALL_COLUMNS[i].replace("fold_", "") for i in FOLD_COLUMNS])
    data_dict[name]["bind_colnames"] = np.asarray([ALL_COLUMNS[i].replace("bind_", "") for i in BIND_COLUMNS])
  return data_dict

#Resample training data
def resample_training_data(tensor_dict, n_resamplings, rand_seed):
  #Resample observed fitness from error distribution
  np.random.seed(rand_seed)
  observed_fitness = np.array(tensor_dict["target"])
  observed_fitness_sd = np.array(tensor_dict["target_sd"])
  observed_fitness_resample = np.array(
    [np.array(
      [observed_fitness[i]+np.random.normal(0, observed_fitness_sd[i]) for i in range(len(observed_fitness))])
    for j in range(n_resamplings)])
  #Save new data
  tensor_dict["target"] = tf.expand_dims(tf.convert_to_tensor(np.ravel(observed_fitness_resample), np.float32), -1)
  tensor_dict["select"] = tf.concat([tensor_dict["select"] for i in range(n_resamplings)], axis = 0)
  tensor_dict["fold"] = tf.sparse.concat(axis = 0, sp_inputs=[tensor_dict["fold"] for i in range(n_resamplings)])
  tensor_dict["bind"] = tf.sparse.concat(axis = 0, sp_inputs=[tensor_dict["bind"] for i in range(n_resamplings)])
  return(tensor_dict)

#Get sequence ID from sequence string
def get_seq_id(sq):
  return ":".join([str(i)+sq[i] for i in range(len(sq))])

#Little function that returns layer index corresponding to layer name
def get_layer_index(model, layername):
  for idx, layer in enumerate(model.layers):
    if layer.name == layername:
      return idx

def shuffle_weights(model, weights=None):
  """Randomly permute the weights in `model`, or the given `weights`.

  This is a fast approximation of re-initializing the weights of a model.

  Assumes weights are distributed independently of the dimensions of the weight tensors
    (i.e., the weights have the same distribution along each dimension).

  :param Model model: Modify the weights of the given model.
  :param list(ndarray) weights: The model's weights will be replaced by a random permutation of these weights.
    If `None`, permute the model's current weights.
  """
  if weights is None:
    weights = model.get_weights()
  weights = [np.random.permutation(w.flat).reshape(w.shape) for w in weights]
  # Faster, but less random: only permutes along the first dimension
  # weights = [np.random.permutation(w) for w in weights]
  model.set_weights(weights)

#Function to create a keras model
def create_model(learn_rate, l1, l2, input_dim_select, input_dim_folding, input_dim_binding, number_additive_traits):
  ### INPUT LAYER
  ########################################
  #Input layer select
  input_layer_select = keras.layers.Input(
    shape = input_dim_select,
    name = "input_select")
  #Split
  input_layer_select_folding = keras.layers.Lambda(lambda x: tf.expand_dims(x[:,0],-1))(input_layer_select)
  input_layer_select_binding = keras.layers.Lambda(lambda x: tf.expand_dims(x[:,1],-1))(input_layer_select)
  #Input layer folding
  input_layer_folding = keras.layers.Input(
    shape = input_dim_folding,
    name = "input_fold", sparse=True)
  #Input layer binding
  input_layer_binding = keras.layers.Input(
    shape = input_dim_binding,
    name = "input_bind", sparse=True)
  ### FOLDING LAYERS
  ########################################
  #Folding additive trait layer
  folding_additive_trait_layer = keras.layers.Dense(
    number_additive_traits,
    # input_dim = input_dim,
    kernel_initializer = 'glorot_normal',
    activation = "linear",
    use_bias = False,
    name = "folding_additivetrait")(input_layer_folding)
    # kernel_regularizer = keras.regularizers.l1_l2(l1 = l1, l2 = l2))(input_layerB)
  #Folding nonlinear layer
  folding_nonlinear_layer = State_prob_folded(trainable=False)(folding_additive_trait_layer)
  #Folding additive layer
  folding_additive_layer = keras.layers.Dense(
    1,
    activation = "linear",
    name = "folding_additive",
    kernel_constraint=Between(0, 1e3)#, bias_constraint=Between(-1, 1)
    )(folding_nonlinear_layer)
  ### BINDING LAYERS
  ########################################
  #Binding additive trait layer
  binding_additive_trait_layer = keras.layers.Dense(
    number_additive_traits,
    # input_dim = input_dim,
    kernel_initializer = 'glorot_normal',
    activation = "linear",
    use_bias = False,
    name = "binding_additivetrait",
    kernel_regularizer = keras.regularizers.l1_l2(l1 = l1, l2 = l2))(input_layer_binding)
  #Binding nonlinear layer
  binding_nonlinear_layer = State_prob_bound(trainable=False)([binding_additive_trait_layer, folding_additive_trait_layer])
  #Binding additive layer
  binding_additive_layer = keras.layers.Dense(
    1,
    activation = "linear",
    name = "binding_additive",
    kernel_constraint=Between(0, 1e3)#, bias_constraint=Between(-1, 1)
    )(binding_nonlinear_layer)
  ### OUTPUT LAYERS
  ########################################
  #Multiplicative layer folding
  multiplicative_layer_folding = keras.layers.Multiply()([folding_additive_layer, input_layer_select_folding])
  #Multiplicative layer binding
  multiplicative_layer_binding = keras.layers.Multiply()([binding_additive_layer, input_layer_select_binding])
  #Sum layer
  output_layer = keras.layers.Add()([multiplicative_layer_folding, multiplicative_layer_binding])
  #Create keras model defining input and output layers
  model = keras.Model(
    inputs = [input_layer_select, input_layer_folding, input_layer_binding], 
    outputs = [output_layer])
  # Compile model
  opt = keras.optimizers.Adam(learning_rate = learn_rate)
  #Compile the model
  model.compile(
    optimizer = opt,
    loss = 'mean_absolute_error')
  return model

#Fit model for gridsearch
def fit_model_grid(param_dict, input_data, n_epochs):
  #Clear session
  keras.backend.clear_session()
  #Summarize results
  print("Grid search using %s" % (param_dict))
  #Set random seeds
  random.seed(random_seed)
  tf.random.set_seed(random_seed)
  #Create model
  model = create_model(
    learn_rate = param_dict['learning_rate'],
    l1=param_dict['l1_regularization_factor'],
    l2=param_dict['l2_regularization_factor'],
    input_dim_select = input_data['train']['select'].shape[1], 
    input_dim_folding = input_data['train']['fold'].shape[1], 
    input_dim_binding = input_data['train']['bind'].shape[1], 
    number_additive_traits = param_dict['number_additive_traits'])
  #Validation data
  validation_data = (
    [input_data['valid']['select'], input_data['valid']['fold'], input_data['valid']['bind']],
    input_data['valid']['target'])
  #Fit the model
  history = model.fit(
    [input_data['train']['select'], input_data['train']['fold'], input_data['train']['bind']],
    input_data['train']['target'],
    validation_data = validation_data,
    epochs = n_epochs,
    batch_size = param_dict['num_samples'],
    shuffle = True,
    verbose = 0,
    use_multiprocessing = True)
  return(history.history["val_loss"][-1])

#######################################################################
## SETUP ##
#######################################################################

#Output model directory
model_directory = os.path.join(output_directory, "whole_model")
#Create output model directory
try:
  os.mkdir(model_directory)
except FileExistsError:
  print("Warning: Output model directory already exists.")

#Output plot directory
plot_directory = os.path.join(output_directory, "plots")
#Create output plot directory
try:
  os.mkdir(plot_directory)
except FileExistsError:
  print("Warning: Output plot directory already exists.")

#Load model data
model_data = load_model_data({
  "train": data_train_file,
  "valid": data_valid_file,
  "obs": data_obs_file})

#Resample training data
if num_resamplings!=0:
  model_data["train"] = resample_training_data(
    tensor_dict = model_data["train"], 
    n_resamplings = num_resamplings, 
    rand_seed = random_seed)

#######################################################################
## TUNE LEARNING RATE, NUMBER OF SAMPLES AND REGULARISATION PARAMS ##
#######################################################################

if len(l1)==1 and len(l2)==1 and len(batch_size)==1 and len(learn_rate)==1:
  #Only parameters
  num_samples = batch_size[0]
  learning_rate = learn_rate[0]
  l1_regularization_factor = l1[0]
  l2_regularization_factor = l2[0]
else:
  #All combinations of tunable parameters
  parameter_grid = [{
  "num_samples":i, 
  "learning_rate":j, 
  "l1_regularization_factor":k, 
  "l2_regularization_factor":l, 
  "number_additive_traits":1} for i in batch_size for j in learn_rate for k in l1 for l in l2]
  #Perform grid search
  grid_results = [fit_model_grid(i, model_data, num_epochs_grid) for i in parameter_grid]
  best_params = parameter_grid[[i for i in range(len(grid_results)) if grid_results[i]==min(grid_results)][0]]
  #Summarize results
  print("Best: %f using %s" % (min(grid_results), best_params))
  #Best parameters
  num_samples = best_params['num_samples']
  learning_rate = best_params['learning_rate']
  l1_regularization_factor = best_params['l1_regularization_factor']
  l2_regularization_factor = best_params['l2_regularization_factor']

#######################################################################
## BUILD FINAL NEURAL NETWORK ##
#######################################################################

#Clear session
keras.backend.clear_session()

#Set random seeds
random.seed(random_seed)
tf.random.set_seed(random_seed)

#Create model
model = create_model(
  learn_rate = learning_rate,
  l1=l1_regularization_factor,
  l2=l2_regularization_factor,
  input_dim_select = model_data['train']['select'].shape[1], 
  input_dim_folding = model_data['train']['fold'].shape[1], 
  input_dim_binding = model_data['train']['bind'].shape[1], 
  number_additive_traits = number_additive_traits)
print(model.summary())

#Validation data
validation_data = (
  [model_data['valid']['select'], model_data['valid']['fold'], model_data['valid']['bind']],
  model_data['valid']['target'])

#Save model weights
original_model_weights = model.get_weights()

#Fit model(s)
for model_count in range(num_models):

  #Shuffle model weights
  shuffle_weights(model, original_model_weights)

  #Callbacks
  model_callbacks = []
  if early_stopping:
    model_callbacks = [
      keras.callbacks.EarlyStopping(monitor='val_loss', patience=num_epochs*0.1), 
      keras.callbacks.ModelCheckpoint(filepath=os.path.join(model_directory, 'my_model_'+str(model_count)), monitor='val_loss', save_best_only=True)]

  #Fit the model
  history = model.fit(
    [model_data['train']['select'], model_data['train']['fold'], model_data['train']['bind']],
    model_data['train']['target'],
    validation_data = validation_data,
    epochs = num_epochs,
    batch_size = num_samples,
    shuffle = True,
    callbacks = model_callbacks,
    verbose = 2,
    use_multiprocessing = True)

  #Save the entire model as a SavedModel
  if early_stopping==False:
    model.save(os.path.join(model_directory, 'my_model_'+str(model_count)))

  #Load model
  custom = {'Between':Between}
  model = load_model(os.path.join(model_directory, 'my_model_'+str(model_count)))

  #Plot model performance per epoch
  my_figure = plt.figure(figsize = (8,8))
  plt.plot(
    np.log(history.history['loss']))
  plt.xlabel('Number of epochs')
  plt.ylabel('Mean Absolute Error (MAE) on testing data')
  my_figure.savefig(os.path.join(plot_directory, "model_performance_perepoch_"+str(model_count)+".pdf"), bbox_inches='tight')

  #######################################################################
  ## SAVE OBSERVATIONS, PREDICTIONS & ADDITIVE TRAIT VALUES ##
  #######################################################################

  #Model predictions on observed variants
  model_predictions = model.predict([model_data['obs']['select'], model_data['obs']['fold'], model_data['obs']['bind']])

  #Index for folding additive trait layer
  layer_idx_folding = get_layer_index(
    model = model,
    layername = "folding_additivetrait")
  #Calculate folding additive trait
  folding_additive_traits_model = keras.Model(
    inputs = model.input,
    outputs = model.layers[layer_idx_folding].output)
  #Convert to data frame
  folding_additive_trait_df = pd.DataFrame(folding_additive_traits_model.predict([model_data['obs']['select'], model_data['obs']['fold'], model_data['obs']['bind']]))
  folding_additive_trait_df.columns = [ "trait " + str(i) for i in range(len(folding_additive_trait_df.columns))]

  #Index for binding additive trait layer
  layer_idx_binding = get_layer_index(
    model = model,
    layername = "binding_additivetrait")
  #Calculate binding additive trait
  binding_additive_traits_model = keras.Model(
    inputs = model.input,
    outputs = model.layers[layer_idx_binding].output)
  #Convert to data frame
  binding_additive_trait_df = pd.DataFrame(binding_additive_traits_model.predict([model_data['obs']['select'], model_data['obs']['fold'], model_data['obs']['bind']]))
  binding_additive_trait_df.columns = [ "trait " + str(i) for i in range(len(binding_additive_trait_df.columns))]

  #Results data frame
  dataframe_to_export = pd.DataFrame({
    "seq" : np.array(model_data['obs']['sequence']).flatten(),
    "observed_fitness" : np.array(model_data['obs']['target']).flatten(),
    "predicted_fitness" : model_predictions.flatten(),
    "additive_trait_folding" : folding_additive_trait_df["trait 0"],
    "additive_trait_binding" : binding_additive_trait_df["trait 0"],
    "training_set" : np.array(model_data['obs']['training_set']).flatten()})
  #Save as csv file
  dataframe_to_export.to_csv(
    os.path.join(output_directory, "predicted_fitness_"+str(model_count)+".txt"),
    sep = "\t",
    index = False)

  #Save model weights
  dataframe_to_export_folding = pd.DataFrame({
    "id" : model_data['obs']['fold_colnames'],
    "folding_coefficient" : [i[0] for i in model.layers[layer_idx_folding].get_weights()[0]]})
  dataframe_to_export_binding = pd.DataFrame({
    "id" : model_data['obs']['bind_colnames'],
    "binding_coefficient" : [i[0] for i in model.layers[layer_idx_binding].get_weights()[0]]})
  #Merge
  dataframe_to_export = dataframe_to_export_folding.merge(dataframe_to_export_binding, left_on='id', right_on='id', how='outer')
  #Save as csv file
  dataframe_to_export.to_csv(
    os.path.join(output_directory, "model_weights_"+str(model_count)+".txt"),
    sep = "\t",
    index = False)

  #Save remaining model parameters (linear layers)
  with open(os.path.join(output_directory, "model_parameters_"+str(model_count)+".txt"), 'w') as f:
    for ml in model.layers:
      if(ml.name in ["folding_additive", "binding_additive"]):
        f.write(ml.name.replace("additive", "linear")+"_kernel\n")
        f.write(str(float(ml.weights[0]))+"\n")
        f.write(ml.name.replace("additive", "linear")+"_bias\n")
        f.write(str(float(ml.weights[1]))+"\n")

