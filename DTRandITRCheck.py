from Bio import SeqIO
import os
import argparse


# Version 3.0


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


def validate_file(f):
    if not os.path.exists(f):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(f))
    return f


# Both matchDTR and match ITR take too much time, somehow merge them together
def matchDTRandITR(a, b):
    """This function checks the percentage to which two DTR and ITR sequences are equal"""
    if len(a) != len(b):
        print("Lengths are not equal. A bug has occurred in the program.")
        exit()
    d = 0
    i = 0
    for y in range(len(a)):
        if a[y] == b[y]:
            d += 1
        if a[len(a) - 1 - y] == b[y]:
            i += 1
    return float(d) / float(len(a)), float(i) / float(len(a))


def maybeSequence(p):
    """This function checks if a percentage is greater than a specified error bound"""
    global errorBound
    return p >= errorBound


def printOutSequences(a, b, l):
    """This is a function that prints out the selected sequences to the user"""
    global sequence
    global IsItDTR
    global DTRBestPercentage
    global ITRBestPercentage
    print('─' * 10)
    if IsItDTR:
        print("Program found the following sequences having a similarity percentage of " + f"{DTRBestPercentage:.2%}"
              + " for being a DTR")
        if maybeSequence(DTRBestPercentage):
            print("This sequence is most likely to be a DTR.")
    else:
        print("Program found the following sequences having a similarity percentage of " + f"{ITRBestPercentage:.2%}"
              + " for being a ITR")
        if maybeSequence(ITRBestPercentage):
            print("This sequence is most likely to be a ITR.")
    print("Sequences have a length of " + str(l) + " character(s)")
    print("The genome file has " + str(len(sequence)) + " characters in total")
    print("Sequence 1: " + sequence[a:a + l])
    print("Starting from index " + str(a) + " to index " + str(a + l - 1))
    print("Sequence 2: " + sequence[b:b + l])
    print("Starting from index " + str(b) + " to index " + str(b + l - 1))


def DTROrITRMatch(a, b, l):
    """This function checks if the sequence similarity is the greatest and updates and returns variables accordingly"""
    global DTRBestPercentage
    global DTRBestPercentage
    global ITRBestPercentage
    global indexFirst
    global indexSecond
    global indexLength
    global sequence
    global IsItDTR
    global DTRFound
    global ITRFound
    DTRSimilarity, ITRSimilarity = matchDTRandITR(sequence[a:a + l], sequence[b:b + l])
    if DTRSimilarity > DTRBestPercentage:
        if DTRSimilarity > ITRBestPercentage:
            IsItDTR = True
        DTRBestPercentage = DTRSimilarity
        indexFirst = a
        indexSecond = b
        indexLength = l
        if maybeSequence(DTRSimilarity):
            DTRFound = True
            printOutSequences(a, b, l)
    if ITRSimilarity > DTRBestPercentage:
        if ITRSimilarity > DTRBestPercentage:
            IsItDTR = False
        ITRBestPercentage = ITRSimilarity
        indexFirst = a
        indexSecond = b
        indexLength = l
        if maybeSequence(ITRSimilarity):
            ITRFound = False
            printOutSequences(a, b, l)


# Getting the path of the genome file
# Currently, this DTRCheck can only access the files of the directory it is currently in <- Rectify in future

# Using argparse to get a path from the argument itself
parser = argparse.ArgumentParser(description="Inputting path of genome file")
sequenceFile = ""
# Add arguments for genome file path
parser.add_argument("-p", "--path", dest="filename", required=False, type=validate_file, help="Genome File Path",
                    metavar="")
args = parser.parse_args()
if args.filename is not None:
    sequenceFile = args.filename
else:
    print("There are two ways to input the genome file")
    print("1: Input a path that leads to the genome file")
    print("2: Search for the genome file within the current directory")
    choice = input("Choice: ")
    if choice == "1":
        sequenceFile = input("Path: ")
    elif choice == "2":
        while True:
            print("Please select a file from the current directory to check.")
            printingFilesInDir()
            index = input("What number file will you select? ")
            ans = gettingSpecificFile(index)
            if ans == "None such file":
                print("Please try again.")
            else:
                sequenceName = ans
                break
    else:
        print("Not a valid choice.")
        exit()

# Loading the sequence from the file
sequence = ""
try:
    sequence = SeqIO.read(sequenceFile, "fasta").seq
    if len(sequence) < 1000:
        print("This is not an acceptable genome file.\nThe length of the sequence is too small")
        exit()
    print("File has been loaded.")
    print('─' * 10)
except ValueError:
    print("Sequence didn't load. This file not a .fasta file.")
    exit()

minimumDTR = 70
try:
    print("What is the minimum length of the DTR/ITR sequence you want to search?")
    print("(If the value is invalid, the value will be set to 70 by default)")
    minimumDTR = int(input("Minimum Length: "))
except ValueError:
    print("Value is considered invalid. Default value (70) was implemented.")
print('─' * 10)

errorBound = 0.95
try:
    print("What is the error bound percentage of the DTR/ITR sequence you want to search?")
    print("(If the value is invalid, the value will be set to 95% by default)")
    errorBound = float(input("Error Bound: ")) / 100
except ValueError:
    print("Value is considered invalid. Default value (95%) was implemented.")
print('─' * 10)

window = 300
try:
    print("What is the window of the ends of the genome file you want to search?")
    print("(If the value is invalid, the value will be set to 300 characters by default)")
    errorBound = float(input("Window: "))
except ValueError:
    print("Value is considered invalid. Default value (300) was implemented.")
print('─' * 10)

print("Calculating...")
indexFirst = -1
indexLength = -1
indexSecond = -1
DTRBestPercentage = -1.0
ITRBestPercentage = -1.0
step = 10
stepCopy = step
DTRFound = False
ITRFound = False
IsItDTR = False

# Split the sequence into two halves of 300 and search for the best matching sequence in each
for MatchLength in range(200, minimumDTR - step, -step):
    print("Checking for DTR/ITR of length: " + str(MatchLength))
    for i in range(window - MatchLength + 1):
        for j in range(len(sequence) - window - 1, len(sequence) - MatchLength + 1):
            DTROrITRMatch(i, j, MatchLength)
            if DTRFound or ITRFound:
                break
        if DTRFound or ITRFound:
            break
    if DTRFound or ITRFound:
        break
    if DTRBestPercentage > ITRBestPercentage:
        IsItDTR = True
    else:
        IsItDTR = False
    print("DTR best percentage found: " + f"{DTRBestPercentage:.2%}")
    print("ITR best percentage found: " + f"{ITRBestPercentage:.2%}")

if not DTRFound and not ITRFound:
    printOutSequences(indexFirst, indexSecond, indexLength)

# Now, expand the indexes
print('─' * 10)
print("Expanding indexes...")


def expandDTRIndex():
    """If DTRPercentage is greater, then this looks on the sides of the sequence to check if any other characters
    match"""
    global indexFirst
    global indexSecond
    global step
    global sequence
    global indexLength
    # Left Index
    while indexFirst > 0 and step > 0:
        if sequence[indexFirst - 1] == sequence[indexSecond - 1]:
            indexFirst -= 1
            indexSecond -= 1
            step -= 1
            indexLength += 1
        else:
            break
    # Right Index
    while indexSecond < len(sequence) - 1 - indexLength and step > 0:
        if sequence[indexFirst + indexLength] == sequence[indexSecond + indexLength]:
            step -= 1
            indexLength += 1
        else:
            break


def expandITRIndex():
    """If ITRPercentage is greater, then this looks on the sides of the sequence to check if any other characters
    match"""
    global indexFirst
    global indexSecond
    global step
    global sequence
    global indexLength
    # Internal Index
    while indexSecond - indexFirst > MatchLength - 2 and step > 0:
        if sequence[indexFirst + indexLength] == sequence[indexSecond - 1]:
            indexSecond -= 1
            step -= 1
            indexLength += 1
        else:
            break
    # External Index
    while indexFirst > 0 and indexSecond < len(sequence) - 1 - indexLength and step > 0:
        if sequence[indexFirst - 1] == sequence[indexSecond + indexLength]:
            indexFirst -= 1
            step -= 1
            indexLength += 1
        else:
            break


if IsItDTR:
    expandDTRIndex()
    DTRBestPercentage, _ = matchDTRandITR(sequence[indexFirst:indexFirst + indexLength],
                                      sequence[indexSecond:indexSecond + indexLength])
else:
    expandITRIndex()
    _, ITRBestPercentage = matchDTRandITR(sequence[indexFirst:indexFirst + indexLength],
                                      sequence[indexSecond:indexSecond + indexLength])
if step == stepCopy:
    print("Sequence was not able to be expanded.")
else:
    print("Sequence was expanded by " + str(stepCopy - step) + " character(s)")
    printOutSequences(indexFirst, indexSecond, indexLength)

# Print the results to an output file
print('─' * 10)
if input("Would you like to output the results to a file? (Say yes if true)\n").lower() == "yes":
    if IsItDTR:
        fileName = "DTROutput.txt"
    else:
        fileName = "ITROutput.txt"
    if input("Is there a specific path you would like to output to? (Say yes if true)\n").lower() == "yes":
        path = input("Input a path here: ")
        fileName = os.path.join(path, fileName)
    else:
        print("You have not selected to output to a specific path")
    output = ""
    if fileName != "":
        output = open(fileName, "a")
    else:
        print("You have not specified a valid path.")
        exit()
    if IsItDTR:
        output.write("Program found the following sequences having a similarity percentage of "
                     + f"{DTRBestPercentage:.2%} for being a DTR\n")
        if maybeSequence(DTRBestPercentage):
            output.write("This sequence is most likely to be a DTR.\n")
    else:
        output.write("Program found the following sequences having a similarity percentage of "
                     + f"{ITRBestPercentage:.2%} for being a ITR\n")
        if maybeSequence(ITRBestPercentage):
            output.write("This sequence is most likely to be a ITR.\n")
    output.write("Sequences have a length of " + str(indexLength) + " character(s)\n")
    output.write("The genome file has " + str(len(sequence)) + " characters in total\n")
    output.write("Sequence 1: " + str(sequence[indexFirst:indexFirst + indexLength]) + "\n")
    output.write("Starting from index " + str(indexFirst) + " to index " + str(indexFirst + indexLength - 1) + "\n")
    output.write("Sequence 2: " + str(sequence[indexSecond:indexSecond + indexLength]) + "\n")
    output.write("Starting from index " + str(indexSecond) + " to index " + str(indexSecond + indexLength - 1) + "\n")
else:
    print("You have not selected to output the results to a file")
exit()

# Possible Improvements:
#   -Have full file traversal
