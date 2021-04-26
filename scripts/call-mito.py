import vcf
import gzip
import pdb
import sys
# raw_variants = vcf.Reader(open('variants/grouped/grouped.ncbiformat.vcf'))
#pdb.set_trace()
#from vcf.model import remake_calldata_tuple
from vcf.model import _Call

infile = sys.argv[1]
outfile = sys.argv[2]

raw_variants = vcf.Reader(filename=infile)

# to add info tag for AF (aka VF?)
#raw_variants.infos['HOB']
# _Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc', 'source', 'version'])

# todo, call _parse_sample_format
raw_variants.formats['AF'] = vcf.parser._Info(id="AF", num=1, type="Float", desc="Allele frequencyfor individual (named inconsistently in literature for site/sample aka VD",source=None, version=None)
# use modified raw_variants as template
vcf_writer = vcf.Writer(open(outfile, 'w'), raw_variants)

 
for record in raw_variants:
#     print(record.INFO)
    # Add AF to site Format

    # this way is closer to how the internals work. I need to hard code for now because designed to be called while parsing samples
    #record.FORMAT = raw_variants._parse_sample_format(":".join([record.FORMAT,"AF"]))

    #record.FORMAT = ":".join([record.FORMAT,"AF"])
    oldFormat = record.FORMAT

    
    record.add_format("AF")


    # copying from _parse_samples method in parser.py
    samp_fmt = raw_variants._parse_sample_format(record.FORMAT)
    
    # _parse_sample_format method in parser.py seems to validate this
    nfields = len(samp_fmt._fields)
    samp_data = []
    
    for sample_call in record:
        
        tempTuple = sample_call.data

        #if(len(tempTuple.AD) > 2):
        #    pdb.set_trace()

        # add AF
        refCount = tempTuple.AD[0]
        varCount = tempTuple.AD[1]
        totalBases = refCount + varCount

        if((refCount + varCount) > 0):
            tempAF = varCount / (totalBases)
        else:
            # should see if there is a valid null for floats
            tempAF = 0.0

        # eventually may add "FT" for genotype filter

        newGT = "0/0"
        # call preliminary genotype
        if(totalBases >= 100):
            if((tempAF > 0.01) and (tempAF < 0.99)):
                newGT = "0/1"
            elif(tempAF >= 0.99):
                newGT = "1/1"
            else:
                newGT = "0/0"
        elif((totalBases < 100) and (totalBases >= 10)):
            # if the total bases is less than 100 we allow homoplasmies
            if(tempAF >= 0.99):
                newGT = "1/1"
            else:
                newGT = "0/0"
        
        # make list of vals in correct order. need to edit GT
        callData = []
        for name in oldFormat.split(":"):
            if(name == "GT"):
                callData.append(newGT)
            else:
            
                tempVal = getattr(sample_call.data, name)
                callData.append(tempVal)
        #callData = [i for i in sample_call.data]
        callData.append(tempAF)
        #if(tempAD >= 0.90):
        #    tempTuple = tempTuple
            #pdb.set_trace()
    

        # process AF and make calls like in parser.py
        # then create call _Call(site, name, samp_fmt(*sampdat))
        
        call = _Call(record, sample_call.sample, samp_fmt(*callData)) 
        samp_data.append(call)
    record.samples = samp_data
    vcf_writer.write_record(record)
        
