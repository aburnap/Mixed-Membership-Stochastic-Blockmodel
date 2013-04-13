import System
from System import Console
from System import Array
from System import IO
import clr
import sys

infernetpackagedir = System.IO.Path.GetDirectoryName(System.IO.FileInfo(__file__).FullName)

#add path to infer.net dlls .................................
clr.AddReferenceToFileAndPath(infernetpackagedir+"\\Infer.Runtime.dll")
clr.AddReferenceToFileAndPath(infernetpackagedir+"\\Infer.Compiler.dll")

#---import all namespaces----------------
import MicrosoftResearch.Infer
import MicrosoftResearch.Infer.Models
import MicrosoftResearch.Infer.Distributions
import MicrosoftResearch.Infer.Collections
import MicrosoftResearch.Infer.Factors
import MicrosoftResearch.Infer.Maths
from MicrosoftResearch.Infer import * 

#---------import all classes and methods from above namespaces-----------
from MicrosoftResearch.Infer.Distributions import *
from MicrosoftResearch.Infer.Models import *
from MicrosoftResearch.Infer.Collections import *
from MicrosoftResearch.Infer.Factors import *
from MicrosoftResearch.Infer.Maths import *

#---import wrapper---
import InferNetWrapper
