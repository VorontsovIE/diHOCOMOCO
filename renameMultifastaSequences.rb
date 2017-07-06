raise 'Specify control name' unless control_name = ARGV[0]
$stdin.readlines.each_slice(2) do |hdr, seq|
  puts ">#{control_name}:#{seq.strip.length}\n#{seq}"
end
