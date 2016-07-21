from utilities import *

#is there a way to add dependencies in dispersion and electr?
#for example, if lj then require cuttype and cutoff, but if none, don't require anything

#can i limit the number of mcf files to the number of species?

#do maxDisp schema

#ad unique to avoid duplication of keywords

#with open('cassandraSchema.json','r') as f:
#	cassandraSchema = json.load(f)
#with open('input.json','r') as f:
#	inputFile = json.load(f)

validateJsonFile("input.json", "cassandraSchema.json")


jsonF = open("input.json",'r')
jsonData = json.load(jsonF)

json2input(jsonData,"input.mcf")

