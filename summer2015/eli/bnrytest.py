def bnrytest(place,num):
    digit1=2**place; digit2=2**(place+1)
    bin=(num%digit2)/digit1
    print bin
bnrytest(7,36)
