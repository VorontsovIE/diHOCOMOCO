ROCPoint = Struct.new(:tpr, :fpr) do # upper boundary for tpr is reported
  def to_s; "<fpr:#{fpr};tpr:#{tpr}>"; end
  def inspect; to_s; end
  def round(round = nil)
    round ? self.class.new(tpr.round(round), fpr.round(round)) : self.class.new(tpr.round, fpr.round)
  end
end

def calculate_auc(roc_points)
  roc_points.each_cons(2).map{|pt_1, pt_2|
    (pt_1.tpr + pt_2.tpr) * (pt_2.fpr - pt_1.fpr) / 2.0
  }.inject(0.0, &:+)
end

def calculate_auc_by_xy(points)
  points.sort_by{|x,y| x }.each_cons(2).map{|pt_1, pt_2|
    (pt_2[0] - pt_1[0]) * (pt_2[1] + pt_1[1]) / 2.0
  }.inject(0.0, &:+)
end

def roc_curve(corrected_pvalues)
  num_positive_sequences = corrected_pvalues.size # all sequences in dataset are positive

  roc_points = []
  roc_points << ROCPoint.new(0.0, 0.0)
  # list pvalues from the most stringent threshold (no one sequence taken as a positive one) to the weakest
  sorted_pvalues = corrected_pvalues.sort
  sorted_pvalues.each_with_index{|corrected_pvalue, index|
    num_positive_sequences_taken = index + 1
    tpr = num_positive_sequences_taken.to_f / num_positive_sequences
    fpr = corrected_pvalue
    roc_points << ROCPoint.new(tpr, fpr)
  }
  roc_points << ROCPoint.new(1.0, 1.0)
  roc_points
end

def precision_recall_curve(corrected_pvalues, vocab_volume)
  num_positive_sequences = corrected_pvalues.size # all sequences in dataset are positive

  roc_points = []
  roc_points << [0.0, 1.0]
  # list pvalues from the most stringent threshold (no one sequence taken as a positive one) to the weakest
  sorted_pvalues = corrected_pvalues.sort
  sorted_pvalues.each_with_index{|corrected_pvalue, index|
    num_positive_sequences_taken = index + 1
    tp = num_positive_sequences_taken.to_f
    recall = tpr = tp / num_positive_sequences
    fpr = corrected_pvalue
    fp = fpr * vocab_volume
    precision = tp / (tp+fp)
    roc_points << [recall, precision]
  }
  roc_points
end
