from Bio import SeqIO
from difflib import SequenceMatcher
import os

# DTR v2.2


def printingFilesInDir():
    """This function lists out all files in the directory"""
    f = 1
    for file in os.listdir():
        print(str(f) + ": " + file)
        f += 1


def gettingSpecificFile(target):
    """This function gets a specified file from the directory"""
    f = 1
    for file in os.listdir():
        if int(f) == int(target):
            return file
        f += 1
    return "None such file"


# Check how SequenceMatcher actually works
# Not counting exactly how similar the two strings are
# Now have made an O(n) function to get the exact result
def match(a, b):
    """This function checks the percentage to which two sequences are equal"""
    if len(a) != len(b):
        print("Lengths are not equal.")
        exit()
    o = 0
    for y in range(len(a)):
        if a[y] == b[y]:
            o += 1
    # m = SequenceMatcher(a=a, b=b).ratio()
    return float(o) / float(len(a))


def maybeDTR(p):
    """This function checks if a percentage is greater than a specified error bound"""
    return p >= 0.95


def printOutSequences(a, b, l, p):
    """This is a function that prints out the selected sequences to the user"""
    global sequence
    print('─' * 10)
    print("Program found the following sequences having a similarity percentage of " + f"{p:.2%}")
    if maybeDTR(p):
        print("This sequence is most likely to be a DTR.")
    print("Sequences have a length of: " + str(l))
    print("The genome file has " + str(len(sequence)) + " characters in total")
    print("Sequence 1: " + sequence[a:a + l])
    print("Starting from index " + str(a) + " to index " + str(a + l - 1))
    print("Sequence 2: " + sequence[b:b + l])
    print("Starting from index " + str(b) + " to index " + str(b + l - 1))


def check(i, j, l):
    """This function checks if the sequence similarity is the greatest and updates and returns variables accordingly"""
    global bestPercentage
    global indexFirst
    global indexSecond
    global indexLength
    global sequence
    similarity = match(sequence[i:i + l], sequence[j:j + l])
    if similarity > bestPercentage:
        bestPercentage = similarity
        indexFirst = i
        indexSecond = j
        indexLength = l
        if maybeDTR(similarity):
            printOutSequences(i, j, l, similarity)
            return True
    return False


# Easy way for a user to get the file from the directory
# Currently, this DTRCheck can only access the files of the directory it is currently in <- Rectify in future
sequenceFile = ""
while True:
    print("Please select a file from the current directory to check.")
    printingFilesInDir()
    index = input("What number file will you select? ")
    ans = gettingSpecificFile(index)
    if ans == "None such file":
        print("Please try again.")
    else:
        sequenceName = ans
        print(ans)
        print("File has been loaded.")
        print('─' * 10)
        break

# Loading the sequence from the file
sequence = ""
try:
    sequence = SeqIO.read(sequenceName, "fasta").seq
    if len(sequence) < 1000:
        print("This is not an acceptable genome file.\nThe length of the sequence is too small")
        exit()
except ValueError:
    print("Sequence didn't load. This file not a .fasta file.")
    exit()


minimumDTR = 70
try:
    print("What is the minimum length of the DTR you want to search?")
    print("(If the value is invalid, the value will be set to 70 by default)")
    minimumDTR = int(input("Minimum Length: "))
except ValueError:
    print("Value is considered invalid.")
    minimumDTR = 70

print("Calculating...")
indexFirst = -1
indexLength = -1
indexSecond = -1
bestPercentage = -1.0
step = 10
stepCopy = step
DTRFound = False

# Split the sequence into two halves of 300 and search for the best matching sequence in each
for DTRLength in range(200, minimumDTR - step, -step):
    print("Checking for DTR of length: " + str(DTRLength))
    for i in range(300 - DTRLength + 1):
        for j in range(len(sequence) - 300 - 1, len(sequence) - DTRLength + 1):
            DTRFound = check(i, j, DTRLength)
            if DTRFound:
                break
        if DTRFound:
            break
    if DTRFound:
        break
    print("Best percentage found: " + f"{bestPercentage:.2%}")
if not DTRFound:
    printOutSequences(indexFirst, indexSecond, indexLength, bestPercentage)

# Now, expand the indexes
print('─' * 10)
print("Expanding indexes...")


def expandFirstIndex():
    """Looks on the left side of the sequence to check if any other characters match"""
    global indexFirst
    global indexSecond
    global step
    global sequence
    global indexLength
    while indexFirst > 0 and step > 0:
        if sequence[indexFirst - 1] == sequence[indexSecond - 1]:
            indexFirst -= 1
            indexSecond -= 1
            step -= 1
            indexLength += 1
        else:
            break


def expandSecondIndex():
    """Looks on the right side of the sequence to check if any other characters match"""
    global indexFirst
    global indexSecond
    global step
    global sequence
    global indexLength
    while indexSecond < len(sequence) - 1 - indexLength and step > 0:
        if sequence[indexFirst + indexLength] == sequence[indexSecond + indexLength]:
            step -= 1
            indexLength += 1
        else:
            break


expandFirstIndex()
expandSecondIndex()
if step == stepCopy:
    print("Sequence was not able to be expanded.")
else:
    print("Sequence was expanded by " + str(stepCopy - step) + " character(s)")
    printOutSequences(indexFirst, indexSecond, indexLength, match(sequence[indexFirst:indexFirst + indexLength], sequence[indexSecond:indexSecond + indexLength]))
exit()

# Possible Improvements:
#   -Might want to consider searching a wider window instead of just 300 characters
#   -Can set up command line arguments using argparse
#   -Have a sequence length cutoff
#       -The length cutoff is usually already in place, just check if the genome has a large enough size
#   -Employ jellyfish.jaro_distance() or Levenshtein distance matching
#   -Have the user decide the error percentage
#   -Have full file traversal
#       -Can also have the user input their own file path
#   -Have an option to write DTR results to an output file
#   -Implementing ITRs
