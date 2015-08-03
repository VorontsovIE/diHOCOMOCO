require 'bioinform'
require 'dipm'
require 'pcm2pwm_converter_fixed'

namespace :collect_and_normalize_data do
  desc 'Convert PCM --> PWM'
  task :convert_pcm_to_pwm do
    mkdir_p 'models/pwm/mono/all'
    FileList['models/pcm/mono/all/*/*.pcm'].each{|fn|
      pcm_to_pwm(fn, fn.pathmap('models/pwm/mono/all/%-1d/%n.pwm'))
    }
    mkdir_p 'models/pwm/di/all'
    FileList['models/pcm/di/all/*/*.dpcm'].each{|fn|
      dipcm_to_dipwm(fn, fn.pathmap('models/pwm/di/all/%-1d/%n.dpwm'))
    }
  end
end
