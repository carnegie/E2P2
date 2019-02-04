from operator import itemgetter

# Classes
class Ensemble(object):
	def __init__(self, sequence_id):
		self.sequence_id = sequence_id
		self.weight = {}
		self.classifiers = {}
		self.predictions = []

	def add_weight(self, weight):
		self.weight = weight

	def add_classifiers(self, classifiers):
		self.classifiers = classifiers

	def add_predictions(self, predictions):
		self.predictions = predictions


class Voting(Ensemble):
	def __init__(self, sequence_id, classifiers):
		super().__init__(sequence_id)
		self.add_classifiers(classifiers)

	def max_weight_voting(self):
		votes = []
		for cname in self.classifiers:
			c = self.classifiers[cname]
			if self.sequence_id in c.predictions:
				temp = []
				# Record classifier's predictions and weights for that ID for voting and for output.
				for ef_class in c.predictions[self.sequence_id]:
					if ef_class in c.weights:
						ef_weight = c.weights[ef_class]
					else:
						ef_weight = "0.000"
					# Record vote as tuple, as this will allow one EF class to have more than one weight, depending
					# on if more than one classifier called it.
					vote = (cname, ef_class, float(ef_weight))
					votes.append(vote)
					# The classifiers attribute will be used in outputting the full results.
					entry = "%s (%s)" % (ef_class, ef_weight)
					temp.append(entry)
				self.classifiers[cname] = temp
		return sorted(votes, key=itemgetter(2), reverse=True)

	def absolute_threshold(self, votes, threshold):
		x = len(votes)
		if x < 1:
			self.predictions.append("NA")
		else:
			# Sort the votes to find the highest weighted vote.
			high_weight = votes[0][2]
			high_class = votes[0][1]
			entry = "%s (%s)" % (high_class, str(high_weight))
			self.predictions.append(entry)

			# Condense the votes so that each EF class appears only once, with its highest weight.
			condensed_votes = {}
			for vote in votes:
				ef_class = vote[1]
				weight = float(vote[2])
				if ef_class in condensed_votes:
					if weight > condensed_votes[ef_class]:
						condensed_votes[ef_class] = weight
				else:
					condensed_votes[ef_class] = weight

			# Iterate through the votes to find all those that pass the threshold.
			t = float(high_weight - threshold)
			# Check to make sure the threshold is not negative.
			if t < 0.0:
				t = 0.0
			for ef_class in condensed_votes:
				if ef_class != high_class:
					weight = condensed_votes[ef_class]
					if float(weight) >= t:
						entry = "%s (%s)" % (ef_class, weight)
						self.predictions.append(entry)
