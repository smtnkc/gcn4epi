import tensorflow as tf
import numpy as np


def masked_softmax_cross_entropy(preds, labels, mask):
    """Softmax cross-entropy loss with masking."""
    loss = tf.nn.softmax_cross_entropy_with_logits(logits=preds, labels=labels)
    mask = tf.cast(mask, dtype=tf.float32)
    mask /= tf.reduce_mean(mask)
    loss *= mask
    return tf.reduce_mean(loss)


def masked_confusion(preds, labels, mask):
    """Confusion matrix with masking."""
    predicted = tf.argmax(preds, 1)
    labels = tf.argmax(labels, 1)
    predicted = tf.cast(predicted, dtype=tf.int32)
    labels = tf.cast(labels, dtype=tf.int32)
    mask = tf.cast(mask, dtype=tf.int32)

    # Count true positives, true negatives, false positives and false negatives.
    tp = tf.count_nonzero(predicted * labels * mask)
    tn = tf.count_nonzero((predicted - 1) * (labels - 1) * mask)
    fp = tf.count_nonzero(predicted * (labels - 1) * mask)
    fn = tf.count_nonzero((predicted - 1) * labels * mask)

    conf_matrix = [[tn, fp], [fn, tp]]
    return conf_matrix


def masked_accuracy(preds, labels, mask):
    """Accuracy with masking."""
    conf = masked_confusion(preds, labels, mask)
    tn = conf[0][0]
    fp = conf[0][1]
    fn = conf[1][0]
    tp = conf[1][1]

    acc = (tp + tn) / (tp + tn + fp + fn)
    return acc


def masked_f1(preds, labels, mask):
    """F1-measure with masking."""
    conf = masked_confusion(preds, labels, mask)
    tn = conf[0][0]
    fp = conf[0][1]
    fn = conf[1][0]
    tp = conf[1][1]

    # Calculate accuracy, precision, recall and F1 score.
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = (2 * precision * recall) / (precision + recall)
    f1 = tf.cast(f1, tf.float32)
    return f1
