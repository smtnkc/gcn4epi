from __future__ import division
from __future__ import print_function

import os
from re import L
import warnings
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
warnings.filterwarnings('ignore')
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import time
import logging
from utils import *
from models import GCN, MLP

# Settings
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('cell_line', 'GM12878', 'Dataset string.')  # 'GM12878', 'K562'
flags.DEFINE_string('cross_cell_line', None, 'Dataset string.')  # 'GM12878', 'K562'
flags.DEFINE_string('model', 'gcn', 'Model string.')  # 'gcn', 'gcn_cheby', 'dense'
flags.DEFINE_float('learning_rate', 0.01, 'Initial learning rate.')
flags.DEFINE_integer('epochs', 200, 'Number of epochs to train.')
flags.DEFINE_integer('hidden1', 16, 'Number of units in hidden layer 1.')
flags.DEFINE_float('dropout', 0.5, 'Dropout rate (1 - keep probability).')
flags.DEFINE_float('weight_decay', 5e-4, 'Weight for L2 loss on embedding matrix.')
flags.DEFINE_integer('early_stopping', 10, 'Tolerance for early stopping (# of epochs).')
flags.DEFINE_integer('max_degree', 3, 'Maximum Chebyshev polynomial degree.')
flags.DEFINE_integer('seed', 42, 'Random seed.')
flags.DEFINE_integer('k_mer', 5, 'K-mer length.')
flags.DEFINE_float('label_rate', 0.2, 'Label rate in [0.2, 0.1, 0.05].')
flags.DEFINE_integer('frag_len', 200, 'Fragment length (base pairs): 0 (disabled) or 200.')

# Set random seed
np.random.seed(FLAGS.seed)
tf.set_random_seed(FLAGS.seed)

# Load data
adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask = load_data(
    FLAGS.cell_line, FLAGS.cross_cell_line, FLAGS.label_rate, FLAGS.k_mer)
n_nodes = features.shape[0]

# Some preprocessing
features = preprocess_features(features)
if FLAGS.model == 'gcn':
    support = [preprocess_adj(adj)]
    num_supports = 1
    model_func = GCN
elif FLAGS.model == 'gcn_cheby':
    support = chebyshev_polynomials(adj, FLAGS.max_degree)
    num_supports = 1 + FLAGS.max_degree
    model_func = GCN
elif FLAGS.model == 'dense':
    support = [preprocess_adj(adj)]  # Not used
    num_supports = 1
    model_func = MLP
else:
    raise ValueError('Invalid argument for model: ' + str(FLAGS.model))

# Define placeholders
placeholders = {
    'support': [tf.sparse_placeholder(tf.float32) for _ in range(num_supports)],
    'features': tf.sparse_placeholder(tf.float32, shape=tf.constant(features[2], dtype=tf.int64)),
    'labels': tf.placeholder(tf.float32, shape=(None, y_train.shape[1])),
    'labels_mask': tf.placeholder(tf.int32),
    'dropout': tf.placeholder_with_default(0., shape=()),
    'num_features_nonzero': tf.placeholder(tf.int32)  # helper variable for sparse dropout
}

# Create model
model = model_func(placeholders, input_dim=features[2][1], logging=True)

# Initialize session
sess = tf.Session()


# Define model evaluation function
def evaluate(features, support, labels, mask, placeholders):
    t_test = time.time()
    feed_dict_val = construct_feed_dict(features, support, labels, mask, placeholders)
    outs_val = sess.run([model.loss, model.accuracy, model.confusion], feed_dict=feed_dict_val)
    return outs_val[0], outs_val[1], outs_val[2], (time.time() - t_test)


# Init variables
sess.run(tf.global_variables_initializer())

cost_val = []

print('Optimization started!')
t_train_start = time.time()
# Train model
for epoch in range(FLAGS.epochs):

    t1 = time.time()
    # Construct feed dictionary
    feed_dict = construct_feed_dict(features, support, y_train, train_mask, placeholders)
    feed_dict.update({placeholders['dropout']: FLAGS.dropout})

    # Training step
    outs = sess.run([model.opt_op, model.loss, model.accuracy], feed_dict=feed_dict)

    # Validation
    cost, acc, confusion, duration = evaluate(features, support, y_val, val_mask, placeholders)

    cost_val.append(cost)

    # Print results
    t2 = time.time()
    print(("Epoch: {:04d}   train_loss = {:.5f}   train_f1 = {:.5f}   val_loss = {:.5f}   " +
           "val_f1 = {:.5f}   time = {:.5f}   conf = {}").format(epoch + 1, outs[1], outs[2], cost, acc, t2-t1, confusion))

    if epoch > FLAGS.early_stopping and cost_val[-1] > np.mean(cost_val[-(FLAGS.early_stopping+1):-1]):
        print("Early stopping...")
        break

t_train_duration = time.time() - t_train_start
print("Optimization Finished in {:.5f}".format(t_train_duration))


# Testing
test_cost, test_acc, test_confusion, test_duration = evaluate(features, support, y_test, test_mask, placeholders)
print("test_loss = {:.5f}   test_F1 = {:.5f}   test_time = {:.5f}".format(test_cost, test_acc, test_duration))

print(test_confusion)

# LOGS

log_dir = "results/{}".format(FLAGS.seed)
if not os.path.isdir(log_dir):
    os.makedirs(log_dir)

if FLAGS.cross_cell_line == None:
    log_name = '{}/{}'.format(log_dir, FLAGS.cell_line)
else:
    log_name = '{}/{}'.format(log_dir, FLAGS.cell_line + '_' + FLAGS.cross_cell_line)

lr = '{:.2f}'.format(FLAGS.label_rate).split('.')[1]
log_file = "{}.txt".format(log_name)
open(log_file, 'w').close() # clear file content
logging.basicConfig(format='%(message)s', filename=log_file,level=logging.DEBUG)
logging.info("Cell-line                  = {}".format(FLAGS.cell_line))
logging.info("Cross Cell-line            = {}".format(FLAGS.cross_cell_line))
logging.info("Fragment length            = {}".format(FLAGS.frag_len))
logging.info("Random seed                = {}".format(FLAGS.seed))
logging.info("Label rate                 = {:.2f}".format(FLAGS.label_rate))
logging.info("K-mer length               = {}".format(FLAGS.k_mer))
logging.info("Total number of nodes      = {}".format(n_nodes))
logging.info("Labeled train nodes (x)    = {}".format(sum(train_mask)))
logging.info("Validation nodes (vx)      = {}".format(sum(val_mask)))
logging.info("Test nodes (tx)            = {}".format(sum(test_mask)))
logging.info("Unlabeled train nodes (ux) = {}".format(n_nodes - (sum(train_mask) + sum(val_mask) + sum(test_mask))))
logging.info("Train epochs               = {}".format(epoch+1))
logging.info("Train F1                   = {:.5f}".format(outs[2]))
logging.info("Train duration             = {:.5f}".format(t_train_duration))
logging.info("Validation F1              = {:.5f}".format(acc))
logging.info("Test F1                    = {:.5f}".format(test_acc))
logging.info("Test cost                  = {:.5f}".format(test_cost))
logging.info("Test duration              = {:.5f}".format(test_duration))
logging.info("Test confusion matrix      = {}".format(test_confusion))
