from utilities import *

#is there a way to add dependencies in dispersion and electr?
#for example, if lj then require cuttype and cutoff, but if none, don't require anything

#can i limit the number of mcf files to the number of species?

#do maxDisp schema

#ad unique to avoid duplication of keywords



inputFile = open("input.json",'r')
inputJson = json.load(inputFile)

schemaF = open("cassandraSchema.json")
schema = json.load(schemaF)

validateJsonFile(inputJson, schema)

json2input(inputJson,"input.mcf")
