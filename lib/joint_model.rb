require_relative 'models'
require_relative 'information_content'

def same_by?(models, &block)
  characteristics = models.map(&block)
  characteristics.all?{|ch| ch == characteristics.first }
end

# Calculate a characteristic for model. Make sure that all models have the same value or raise.
def take_and_check_consistency(models, &block)
  raise 'Can\'t take characteristic for empty model list'  if models.empty?
  if same_by?(models, &block)
    block.call(models.first)
  else
    raise 'Inconsistent characteristics for joint models #{models.inspect}'
  end
end
