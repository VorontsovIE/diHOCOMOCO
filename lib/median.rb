def median(arr)
  raise 'Median of an empty array is undefined'  if arr.empty?
  sorted_arr = arr.sort
  if sorted_arr.size.odd?
    center_index = sorted_arr.size / 2
    sorted_arr[center_index]
  else
    center_index_ceiled = sorted_arr.size / 2
    (sorted_arr[center_index_ceiled] + sorted_arr[center_index_ceiled - 1]) / 2.0
  end
end

def mean(values)
  raise 'Mean of an empty array is undefined'  if values.empty?
  values.inject(0.0, &:+) / values.size
end
