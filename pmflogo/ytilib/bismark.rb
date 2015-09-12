#!/usr/bin/ruby

module Ytilib
require "rexml/document"
include REXML

class Bismark < Document
  
  def initialize(source = nil, add_dtd = false)
    dtd = add_dtd ? "<!DOCTYPE smallbismark SYSTEM 'smallbismark.dtd'>#{$/}" : ""
    source == nil ? super("<?xml version='1.0' encoding='UTF-8'?>#{$/}#{dtd}") : super(source)
    super(IO.read(source)) if source != nil && root == nil
    if source == nil
      self.add_element("smallbismark")
      # xmlns breaks XPath for a REXML library under Linux, strange, indeed
      # self.add_element("smallbismark", {"xmlns" => "http://bioinform.imb.ac.ru/smallBiSMark/smallbismark.dtd"})
      self.root.add_element("comment", {"name" => "WARNING"}).add_text("This is a draft version of small-BiSMark. Specification is the subject to change!")
    end
  end
  
  def getXML
    beautify
    s = ""; write(s, 1, true)
    s.rstrip!
    return s
  end
  alias get_xml getXML
  
  def get_pm(xpath)
    pwmnode = self.elements[xpath]
    pm = PM.new_pm(pwmnode.attribute("length").value.to_i)
    toi = pwmnode.name == "PCM"
    pwmnode.elements.each("pm-column") { |c|
      position = c.attribute("position").value.to_i - 1
      weights = [c.elements["a"].get_text.value.strip.to_f, 
                c.elements["c"].get_text.value.strip.to_f, 
                c.elements["g"].get_text.value.strip.to_f, 
                c.elements["t"].get_text.value.strip.to_f]
      weights.collect { |w| w.to_i } if toi
      pm['A'][position], pm['C'][position], pm['G'][position], pm['T'][position] = weights[0], weights[1], weights[2], weights[3]
    }
    return pm
  end
  
private
  CONTAIN_NO_TEXT = { 
                      "segment" => :vasya_shmyak,
                      "group" => :vasya_shmyak,
                      "smallbismark" => :vasya_shmyak,
                      "motif" => :vasya_shmyak,
                      "PWM" => :vasya_shmyak,
                      "PCM" => :vasya_shmyak,
                      "PPM" => :vasya_shmyak,
                      "source" => :vasya_shmyak,
                      "factor"  => :vasya_shmyak,
                      "pm-column" => :vasya_shmyak,
                      "word-list" => :vasya_shmyak}
  
  def beautify(node = self)
    if node == self 
      self.delete_if { |e| e.is_a?(Text) }
      self.each { |e| beautify(e) }
    else
      node.delete_if { |e| e.is_a?(Text) } if node.respond_to?(:delete_if) && Bismark::CONTAIN_NO_TEXT.has_key?(node.name)
      node.each { |e| beautify(e) } if node.respond_to?(:each)
    end 
  end

end

end