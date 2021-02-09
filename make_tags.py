
text = []
for l in open('reconstructions.tsv','r'):
    text.append(l.strip().split('\t'))

f = open('probability-tags.tex','w')
for l in text[1:]:
    feat=l[1].replace(' ','.')
    #prob='{:.3f}'.format(float(l[2]))
    prob=l[2][:5]
    print("%<*{}>\n{}%</{}>".format(feat,prob,feat),file=f)


f.close()