minlength = 98
maxlength = 100
targetlength = "98"
outfile_R1 = open("F:\\Code\\python\\filterFastq_Piardend\\1_filter.fastq", 'w')

outfile_R2 = open("F:\\Code\\python\\filterFastq_Piardend\\2_filter.fastq", 'w')


with open("F:\\Code\\python\\filterFastq_Piardend\\1.fastq", 'r') as infile_R1, open("F:\\Code\\python\\filterFastq_Piardend\\2.fastq", 'r') as infile_R2 :
    for line_R1, line_R2  in zip(infile_R1, infile_R2):
        try:
            headinfo_R1 = line_R1
            headinfo_R2 = line_R2
            
            sequence_R1 = infile_R1.__next__().strip()
            sequence_R2 = infile_R2.__next__().strip()

            comment_R1 = infile_R1.__next__()
            comment_R2 = infile_R2.__next__()

            quality_R1 = infile_R1.__next__()
            quality_R2 = infile_R2.__next__()
        except:
            break


        seqlength_R1 = len(sequence_R1)
        seqlength_R2 = len(sequence_R2)

        if ( int(minlength) <= seqlength_R1 <= int(maxlength) and int(minlength) <= seqlength_R2 <= int(maxlength) ):
                headinfo_R1split = headinfo_R1.split(" ")
                headinfo_R1af = headinfo_R1split[0] + " " + headinfo_R1split[1] + " " + "length=" + targetlength + "\n"
                sequence_R1af = sequence_R1[0:int(targetlength)] + "\n"
                comment_R1split = comment_R1.split(" ")
                comment_R1af = comment_R1split[0] + " " + comment_R1split[1] + " " + "length=" + targetlength + "\n"
                quality_R1af = quality_R1[0:int(targetlength)] + "\n"

                headinfo_R2split = headinfo_R2.split(" ")
                headinfo_R2af = headinfo_R2split[0] + " " + headinfo_R2split[1] + " " + "length=" + targetlength+"\n"
                sequence_R2af = sequence_R2[0:int(targetlength)] + "\n"
                comment_R2split = comment_R2.split(" ")
                comment_R2af = comment_R2split[0] + " " + comment_R2split[1] + " " + "length=" + targetlength+"\n"
                quality_R2af = quality_R2[0:int(targetlength)] + "\n"

                outfile_R1.write(headinfo_R1af)
                outfile_R1.write(sequence_R1af)
                outfile_R1.write(comment_R1af)
                outfile_R1.write(quality_R1af)

                outfile_R2.write(headinfo_R2af)
                outfile_R2.write(sequence_R2af)
                outfile_R2.write(comment_R2af)
                outfile_R2.write(quality_R2af)
            
        else:
            pass





