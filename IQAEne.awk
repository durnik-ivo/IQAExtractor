#############################################################################
#
# IQAEne
#
# Script for extracting interaction energies from AIMALL calculation
#
# Tutorial:
#       1) Fill in user edit section
#       2) Run "awk -f IQAEne.awk" in the terminal
#       3) For output into file run as "awk -f IQAExtractor.awk > file.txt"
#       4) Carefully check results, especially your definition of fragments and sum files
#
# Author: Ivo Durnik
# Contact: durnik@mail.muni.cz
#
#############################################################################
BEGIN {

#############################################################################
# START OF USER EDIT SECTION


# Specify Sum files
#
# AB
ABSumFile = "examples/at.sum"

# A(Complex,Free) / B(Complex,Free)
AFreeSumFile = "examples/t.sum"
BFreeSumFile = "examples/a.sum"

# A(Opt,Free) / B(Opt,Free)
AOptSumFile = "examples/t_opt.sum"
BOptSumFile = "examples/a_opt.sum"

# Specify atom numbers inside ABSumFile
#
# Syntax: "1,2,3" OR "1-3" OR "1-2,3"
fragmentA = "1-15"
fragmentB = "16-30"

# Produce xyz files to check fragment selection? (yes = 1 / no = 0)
#
# If control=1, script will write xyz files extracted from AIMALL sum files 
# (AB.xyz, A.xyz, AFree.xyz, AOpt.xyz, B.xyz, BFree.xyz, BOpt.xyz)
# Complex is splited according to user selection
control = 1;

# Change reference state, default "all" to consider all atoms in Free and Opt sum files, 
# don't change unless you know what you are doing
#
# Syntax: "all" OR "1,2,3" OR "1-3" OR "1-2,3"
refA = "all"
refB = "all"

# END OF USER EDIT SECTION
#############################################################################

#############################################################################
# PREPARE ANALYSIS
#############################################################################

# a.u. -> Angstrom
au = 0.529177249;
# a.u. -> kcal/mol
kcal = 627.509608;

# Helper variables
currentFile = "";
i = 0;
fileN = 1;
myNR = 1;

N = 1;

NA = 1;
NB = 1;

NAFree = 1;
NBFree = 1;

NACountFree = 1;
NBCountFree = 1;

NAOpt = 1;
NBOpt = 1;

NACountOpt = 1;
NBCountOpt = 1;

pos = 0;
pos2 = 0;
bool1 = 0;
bool2 = 0;
bool3 = 0;
bool4 = 0;

# Tell awk which files to process (order matters)
ARGV[1] = ABSumFile;
ARGV[2] = AFreeSumFile;
ARGV[3] = BFreeSumFile;
ARGV[4] = AOptSumFile;
ARGV[5] = BOptSumFile;
ARGC = 6;


# Split atom numbers defining molecular fragments into arrays
NsetA = split(fragmentA, setA, ",");
for(i in setA) {
    NSubSet = split(setA[i], subSet, "-");
    if(NSubSet > 1) {
        for(j = subSet[1]; j <= subSet[2]; j++) {
            if(j == subSet[1]){
                delete setA[i];
                NsetA++
            }
            setA[NsetA] = j;
            NsetA++;
        } 
    }
}

NsetB = split(fragmentB, setB, ",");
for(i in setB) {
   NSubSet = split(setB[i], subSet, "-");
    if(NSubSet > 1) {
        for(j = subSet[1]; j <= subSet[2]; j++) {
            if(j == subSet[1]){
                delete setB[i];
                NsetB++
            }
            setB[NsetB] = j;
            NsetB++;
        } 
    }
}

# Split atom numbers for reference states
if(refA != "all") {
    NsetRefA = split(refA, setRefA, ",");
    for(i in setRefA) {
        NSubSet = split(setRefA[i], subSet, "-");
        if(NSubSet > 1) {
            for(j = subSet[1]; j <= subSet[2]; j++) {
                if(j == subSet[1]) {
                    delete setRefA[i];
                    NsetRefA++
                }
                setRefA[NsetRefA] = j;
                NsetRefA++;
            } 
        }
    }
}

if(refB != "all") {
    NsetRefB = split(refB, setRefB, ",");
    for(i in setRefB) {
       NSubSet = split(setRefB[i], subSet, "-");
        if(NSubSet > 1) {
            for(j = subSet[1]; j <= subSet[2]; j++) {
                if(j == subSet[1]) {
                    delete setRefB[i];
                    NsetRefB++
                }
                setRefB[NsetRefB] = j;
                NsetRefB++;
            } 
        }
    }
}

# Check if fragments overlap
for(i in setA) {
    for(j in setB) {
        if( setA[i]==setB[j] ) {
            printf("Fragment selection overlaps!\n");
            printf("A: %s\n", fragmentA);
            printf("B: %s\n", fragmentB);
            exitSwitch = 1;
            exit 1;
        }
    }
}

# Check if control set properly
if( control!=1 && control!=0 ) {
    printf("Control output not set properly!\n");
    printf("1, or 0 expected, but was \"%s\" \n", control);
    exitSwitch = 1;
    exit 1;
}

}

#############################################################################
# COLLECT DATA
#############################################################################
{

#
# Process AB Complex Sum File
#
if(fileN == 1) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:") {
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        AB[N] = $1;
        ABx[N] = $3;
        ABy[N] = $4;
        ABz[N] = $5;
        for(atom in setA) {
            if(setA[atom] == N) {
                A[NA] = $1;
                Ax[NA] = $3;
                Ay[NA] = $4;
                Az[NA] = $5;
                NA++;
                break;
            }
        }
        for(atom in setB) {
            if(setB[atom] == N) {
                B[NB] = $1;
                Bx[NB] = $3;
                By[NB] = $4;
                Bz[NB] = $5;
                NB++;
                break;
            }
        }
        N++;
    }
    if(NF == 0) {
        bool1=0;
    }
    
    # Get Self Energy
    if($0 == "IQA Intraatomic (\"Self\") Energy Components:"){
        bool3=1;
        pos=myNR;
    }
    if( bool3 == 1 && myNR > pos+12){
        posA = 0;
        atomA = -1;
        for( atom in A) {
            if(A[atom] == $1){
                posA = 1;
                atomA = atom;
                break;
            }
        }

        posB = 0;
        atomB = -1;
        for(atom in B){
            if(B[atom] == $1){
                posB = 1;
                atomB = atom;
                break;
            }
        }
        if (posA == 1) {
            ASelfT[atomA] = $3;
            ASelfVne[atomA] = $4;
            ASelfVeeC[atomA] = $6;
            ASelfVeeX[atomA] = $7;
            ASelfVel[atomA] = $4+$6;
        }
        if (posB == 1) {
            BSelfT[atomB] = $3;
            BSelfVne[atomB] = $4;
            BSelfVeeC[atomB] = $6;
            BSelfVeeX[atomB] = $7;
            BSelfVel[atomB] = $4+$6;
        }
    }
    if( bool3 == 1 && myNR == pos+12+N-1 ) {
        bool3 = 0;
    }

    # Get Interfragment Interaction
    if($0 == "IQA Diatomic \"Interaction\" Energy Components:"){
        bool2 = 1;
        pos = myNR;

        # Initialize variables for first cycle
        for(i=1; i<NA; i++) {
            ABVne[i] = 0.0;
            ABVen[i] = 0.0;
            ABVnn[i] = 0.0;
            ABVeeC[i] = 0.0;
            ABVeeX[i] = 0.0;
            ABVel[i] = 0.0;
        }
    }
    if( bool2 == 1 && NF != 0 && myNR > pos+13){

        posA = 0;
        atomA = -1;
        for( atom in A) {
            if(A[atom] == $1){
                posA = 1;
                atomA = atom;
            } else if(A[atom] == $2) {
                posA = 2;
                atomA = atom;
            }
        }

        posB = 0;
        atomB = -1;
        for(atom in B){
            if(B[atom] == $1){
                posB = 1;
                atomB = atom;
            } else if(B[atom] == $2){
                posB = 2;
                atomB = atom;
            }
        }
        if( (posA==1 && posB==2) || (posB==1 && posA==2) ) {
            ABVne[atomA] += $4;
            ABVen[atomA] += $5;
            ABVnn[atomA] += $7;
            ABVeeC[atomA] += $8;
            ABVeeX[atomA] += $9;
            ABVel[atomA] += $4+$5+$7+$8;
        }
    }
    if( NF == 0 ) {
        bool2 = 0;
    }
    
    #Get Intrafragment interaction
    if( $0 == "IQA Atomic Contributions to Diatomic \"Interaction\" Energy Components:" ){
        bool4 = 1;
        pos = myNR;
        sectionCount = 1;

        # Initialize for first cycle
        for(i=1; i<NA; i++) {       
            AVne[i] = 0.0;
            AVen[i] = 0.0;
            AVnn[i] = 0.0;
            AVeeC[i] = 0.0;
            AVeeX[i] = 0.0;
            AVel[i] = 0.0;
        }
        for(i=1; i<NB; i++) {
            BVne[i] = 0.0;
            BVen[i] = 0.0;
            BVnn[i] = 0.0;
            BVeeC[i] = 0.0;
            BVeeX[i] = 0.0;
            BVel[i] = 0.0;
        }
    }
    if( bool4 == 1 && myNR == pos+30 ) {
        read = 1;
    }
    if( bool4 == 1 && read == 1 && NF == 9 && sectionCount <= N){
        posA1 = 0;
        posA2 = 0;
        atomA1 = -1;
        atomA2 = -1;
        for( atom in A) {
            if(A[atom] == $1){
                posA1 = 1;
                atomA1 = atom;
            }
            if(A[atom] == $2) {
                posA2 = 1;
                atomA2 = atom;
            }
        }

        posB1 = 0;
        posB2 = 0;
        atomB1 = -1;
        atomB2 = -1;
        for(atom in B){
            if(B[atom] == $1){
                posB1 = 1;
                atomB1 = atom;
            } 
            if(B[atom] == $2){
                posB2 = 1;
                atomB2 = atom;
            }
        }
        if( posA1 == 1 && posA2 == 1 ) {
            AVne[atomA1] += $4*2;
            AVen[atomA1] += $5*2;
            AVnn[atomA1] += $7*2;
            AVeeC[atomA1] += $8*2;
            AVeeX[atomA1] += $9*2;
            AVel[atomA1] += ($4+$5+$7+$8)*2;
        }
        if( posB1 == 1 && posB2 == 1 ) {
            BVne[atomB1] += $4*2;
            BVen[atomB1] += $5*2;
            BVnn[atomB1] += $7*2;
            BVeeC[atomB1] += $8*2;
            BVeeX[atomB1] += $9*2;
            BVel[atomB1] += ($4+$5+$7+$8)*2;
        }      
    }
    if( bool4 == 1 && read == 1 && NF == 1 ) {
        read = 0;
        pos2 = myNR;
        sectionCount++;
    }
    if( bool4 == 1 && myNR == pos2+6 ) {
        read = 1;
    }
    if( sectionCount == N ) {
        bool4 = 0;
    }

#
# Process A Free Sum File
#
} else if (fileN == 2) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        if(refA == "all") {
            AFree[NAFree] = $1;
            AFreex[NAFree] = $3;
            AFreey[NAFree] = $4;
            AFreez[NAFree] = $5;
            NAFree++;
        } else {
            for(atom in setRefA) {
                if(setRefA[atom] == NACountFree) {
                    AFree[NAFree] = $1;
                    AFreex[NAFree] = $3;
                    AFreey[NAFree] = $4;
                    AFreez[NAFree] = $5;
                    NAFree++;
                    break;
                }
            }
        }
        NACountFree++;
    }

    if(NF == 0 ) {
        bool1 = 0;
    }

    # Get Self Energy
    if($0 == "IQA Intraatomic (\"Self\") Energy Components:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+12){
        atomA = -1
        for(atom in AFree){
            if(AFree[atom] == $1){
                AFreeSelfT[atom] = $3;
                AFreeSelfVne[atom] = $4;
                AFreeSelfVeeC[atom] = $6;
                AFreeSelfVeeX[atom] = $7;
                AFreeSelfVel[atom] = $4+$6;
            }
        }
        #AFreeSelfT[atomA] = $3;
        #AFreeSelfVne[atomA] = $4;
        #AFreeSelfVeeC[atomA] = $6;
        #AFreeSelfVeeX[atomA] = $7;
        #AFreeSelfVel[atomA] = $4+$6;
    }
    if( bool3 == 1 && myNR == pos+12+NAFree-1 ) {
        bool3 = 0;
    }

    #Get Intrafragment interaction
    if( $0 == "IQA Atomic Contributions to Diatomic \"Interaction\" Energy Components:" ){
        bool4 = 1;
        pos = myNR;
        sectionCount = 1;

        # Initialize for first cycle
        for(i=1; i<NA; i++) {       
            AFreeVne[i] = 0.0;
            AFreeVen[i] = 0.0;
            AFreeVnn[i] = 0.0;
            AFreeVeeC[i] = 0.0;
            AFreeVeeX[i] = 0.0;
            AFreeVel[i] = 0.0;
        }
    }
    if( bool4 == 1 && myNR == pos+30 ) {
        read = 1;
    }
    if( bool4 == 1 && read == 1 && NF == 9 && sectionCount <= N){
        posA1 = 0;
        posA2 = 0;
        atomA1 = -1;
        atomA2 = -1;
        atomA = -1;
        for( atom in AFree) {
            if(AFree[atom] == $1){
                posA1 = 1;
                atomA1 = atom;
            }
            if(AFree[atom] == $2) {
                posA2 = 2;
                atomA2 = atom;
            }
        }
        if( posA1 == 1 && posA2 == 2 ) {
            AFreeVne[atomA1] += $4*2;
            AFreeVen[atomA1] += $5*2;
            AFreeVnn[atomA1] += $7*2;
            AFreeVeeC[atomA1] += $8*2;
            AFreeVeeX[atomA1] += $9*2;
            AFreeVel[atomA1] += ($4+$5+$7+$8)*2;
        }     
    }
    if( bool4 == 1 && read == 1 && NF == 1 ) {
        read = 0;
        pos2 = myNR;
        sectionCount++;
    }
    if( bool4 == 1 && myNR == pos2+6 ) {
        read = 1;
    }
    if( sectionCount == NA ) {
        bool4 = 0;
    }
#
# Process B Free Sum File
#
} else if (fileN == 3) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        if(refB == "all") {
            BFree[NBFree] = $1;
            BFreex[NBFree] = $3;
            BFreey[NBFree] = $4;
            BFreez[NBFree] = $5;
            NBFree++;
        } else {
            for(atom in setRefB) {
                if(setRefB[atom] == NBCountFree) {
                    BFree[NBFree] = $1;
                    BFreex[NBFree] = $3;
                    BFreey[NBFree] = $4;
                    BFreez[NBFree] = $5;
                    NBFree++;
                    break;
                }
            }
        }
        NBCountFree++;
    }
    if(NF == 0) {
        bool1=0;
    }

    # Get Self Energy
    if($0 == "IQA Intraatomic (\"Self\") Energy Components:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+12){
        for(atom in BFree){
            if(BFree[atom] == $1){
                atomB = atom;
            }
        }

        BFreeSelfT[atomB] = $3;
        BFreeSelfVne[atomB] = $4;
        BFreeSelfVeeC[atomB] = $6;
        BFreeSelfVeeX[atomB] = $7;
        BFreeSelfVel[atomB] = $4+$6;
    }
    if( bool3 == 1 && myNR == pos+12+NBFree-1 ) {
        bool3 = 0;
    }

    #Get Intrafragment interaction
    if( $0 == "IQA Atomic Contributions to Diatomic \"Interaction\" Energy Components:" ){
        bool4 = 1;
        pos = myNR;
        sectionCount=1;

        # Initialize for first cycle
        for(i=1; i<NBFree; i++) {   
            BFreeVne[i] = 0.0;
            BFreeVen[i] = 0.0;
            BFreeVnn[i] = 0.0;
            BFreeVeeC[i] = 0.0;
            BFreeVeeX[i] = 0.0;
            BFreeVel[i] = 0.0;
        }

    }

    if( bool4==1 && myNR==pos+30 ) {
        read = 1;
    }
    if( bool4 == 1 && read == 1 && NF == 9 && sectionCount <= N){
        posB1 = 0;
        posB2 = 0;
        atomB1 = -1;
        atomB2 = -1;
        for( atom in BFree) {
            if(BFree[atom] == $1){
                posB1 = 1;
                atomB1 = atom;
            }
            if(BFree[atom] == $2) {
                posB2 = 2;
                atomB2 = atom;
            }
        }
        if( posB1 == 1 && posB2 == 2 ) {
            BFreeVne[atomB1] += $4*2;
            BFreeVen[atomB1] += $5*2;
            BFreeVnn[atomB1] += $7*2;
            BFreeVeeC[atomB1] += $8*2;
            BFreeVeeX[atomB1] += $9*2;
            BFreeVel[atomB1] += ($4+$5+$7+$8)*2;
        }     
    }
    if( bool4 == 1 && read == 1 && NF == 1 ) {
        read = 0;
        pos2 = myNR;
        sectionCount++;
    }
    if( bool4 == 1 && myNR == pos2+6 ) {
        read = 1;
    }
    if( sectionCount == NBFree ) {
        bool4 = 0;
    }
#
# Process A(Opt) Sum File
# 
} else if (fileN == 4) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }
    if(bool1 == 1 && NF != 0 && myNR > pos+3) {
        if(refA == "all") {
            AOpt[NAOpt] = $1;
            AOptx[NAOpt] = $3;
            AOpty[NAOpt] = $4;
            AOptz[NAOpt] = $5;
            NAOpt++;
        } else {
            for(atom in setRefA) {
                if(setRefA[atom] == NACountOpt) {
                    AOpt[NAOpt] = $1;
                    AOptx[NAOpt] = $3;
                    AOpty[NAOpt] = $4;
                    AOptz[NAOpt] = $5;
                    NAOpt++;
                    break;
                }
            }
        }
        NACountOpt++;
    }
    if(NF == 0 ) {
        bool1 = 0;
    }

    # Get Self Energy
    if($0 == "IQA Intraatomic (\"Self\") Energy Components:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+12){
        for(atom in AOpt){
            if(AOpt[atom] == $1){
                atomA = atom;
            }
        }
        AOptSelfT[atomA] = $3;
        AOptSelfVne[atomA] = $4;
        AOptSelfVeeC[atomA] = $6;
        AOptSelfVeeX[atomA] = $7;
        AOptSelfVel[atomA] = $4+$6;
    }
    if( bool3 == 1 && myNR == pos+12+NAOpt-1 ) {
        bool3 = 0;
    }

    #Get Intrafragment interaction
    if( $0 == "IQA Atomic Contributions to Diatomic \"Interaction\" Energy Components:" ){
        bool4 = 1;
        pos = myNR;
        sectionCount = 1;

        # Initialize for first cycle
        for(i=1; i<NA; i++) {       
            AOptVne[i] = 0.0;
            AOptVen[i] = 0.0;
            AOptVnn[i] = 0.0;
            AOptVeeC[i] = 0.0;
            AOptVeeX[i] = 0.0;
            AOptVel[i] = 0.0;
        }
    }
    if( bool4 == 1 && myNR == pos+30 ) {
        read = 1;
    }
    if( bool4 == 1 && read == 1 && NF == 9 && sectionCount <= N){
        posA1 = 0;
        posA2 = 0;
        atomA1 = -1;
        atomA2 = -1;
        atomA = -1;
        for( atom in AOpt) {
            if(AOpt[atom] == $1){
                posA1 = 1;
                atomA1 = atom;
            }
            if(AOpt[atom] == $2) {
                posA2 = 2;
                atomA2 = atom;
            }
        }
        if( posA1 == 1 && posA2 == 2 ) {
            AOptVne[atomA1] += $4*2;
            AOptVen[atomA1] += $5*2;
            AOptVnn[atomA1] += $7*2;
            AOptVeeC[atomA1] += $8*2;
            AOptVeeX[atomA1] += $9*2;
            AOptVel[atomA1] += ($4+$5+$7+$8)*2;
        }     
    }
    if( bool4 == 1 && read == 1 && NF == 1 ) {
        read = 0;
        pos2 = myNR;
        sectionCount++;
    }
    if( bool4 == 1 && myNR == pos2+6 ) {
        read = 1;
    }
    if( sectionCount == NA ) {
        bool4 = 0;
    }
#
# Process B(Opt) Sum File
# 
} else if (fileN == 5) {
    # Get Atom Info
    if($0 == "Nuclear Charges and Cartesian Coordinates:"){
        bool1 = 1;
        pos = myNR;
    }

    if(bool1 == 1 && NF != 0 && myNR > pos+3){
        if(refB == "all") {
            BOpt[NBOpt] = $1;
            BOptx[NBOpt] = $3;
            BOpty[NBOpt] = $4;
            BOptz[NBOpt] = $5;
            NBOpt++;
        } else {
            for(atom in setRefB) {
                if(setRefB[atom] == NBCountOpt) {
                    BOpt[NBOpt] = $1;
                    BOptx[NBOpt] = $3;
                    BOpty[NBOpt] = $4;
                    BOptz[NBOpt] = $5;
                    NBOpt++;
                    break;
                }
            }
        }
        NBCountOpt++;
    }

    if(NF == 0) {
        bool1=0;
    }

    # Get Self Energy
    if($0 == "IQA Intraatomic (\"Self\") Energy Components:"){
        bool3 = 1;
        pos = myNR;
    }
    if( bool3 == 1 && myNR > pos+12){
        for(atom in BOpt){
            if(BOpt[atom] == $1){
                atomB = atom;
            }
        }

        BOptSelfT[atomB] = $3;
        BOptSelfVne[atomB] = $4;
        BOptSelfVeeC[atomB] = $6;
        BOptSelfVeeX[atomB] = $7;
        BOptSelfVel[atomB] = $4+$6;
    }
    if( bool3 == 1 && myNR == pos+12+NBOpt-1 ) {
        bool3 = 0;
    }

    #Get Intrafragment interaction
    if( $0 == "IQA Atomic Contributions to Diatomic \"Interaction\" Energy Components:" ){
        bool4 = 1;
        pos = myNR;
        sectionCount=1;

        # Initialize for first cycle
        for(i=1; i<NBOpt; i++) {   
            BOptVne[i] = 0.0;
            BOptVen[i] = 0.0;
            BOptVnn[i] = 0.0;
            BOptVeeC[i] = 0.0;
            BOptVeeX[i] = 0.0;
            BOptVel[i] = 0.0;
        }

    }

    if( bool4==1 && myNR==pos+30 ) {
        read = 1;
    }
    if( bool4 == 1 && read == 1 && NF == 9 && sectionCount <= N){
        posB1 = 0;
        posB2 = 0;
        atomB1 = -1;
        atomB2 = -1;
        for( atom in BOpt) {
            if(BOpt[atom] == $1){
                posB1 = 1;
                atomB1 = atom;
            }
            if(BOpt[atom] == $2) {
                posB2 = 2;
                atomB2 = atom;
            }
        }
        if( posB1 == 1 && posB2 == 2 ) {
            BOptVne[atomB1] += $4*2;
            BOptVen[atomB1] += $5*2;
            BOptVnn[atomB1] += $7*2;
            BOptVeeC[atomB1] += $8*2;
            BOptVeeX[atomB1] += $9*2;
            BOptVel[atomB1] += ($4+$5+$7+$8)*2;
        }     
    }
    if( bool4 == 1 && read == 1 && NF == 1 ) {
        read = 0;
        pos2 = myNR;
        sectionCount++;
    }
    if( bool4 == 1 && myNR == pos2+6 ) {
        read = 1;
    }
    if( sectionCount == NBOpt ) {
        bool4 = 0;
    }
#
# In case of error
#
} else {
    print "Error while processing files!"
    exitSwitch = 1;
    exit 1;
}
myNR++;
}

# Change file
ENDFILE {
    fileN++;
    myNR = 1;
}

#############################################################################
# CALCULATE AND PRINT RESULTS
#############################################################################
END {
# Check for error during data colection
if( exitSwitch == 1) {
    exit 1;
}

# Initialize helper variables
#
# AB
totABVne = 0.0;
totABVen = 0.0;
totABVnn = 0.0;
totABVeeC = 0.0;
totABVeeX = 0.0;
totABVel = 0.0;
totABDI = 0.0;

# A(Complex, B)
totAVne = 0.0;
totAVen = 0.0;
totAVnn = 0.0;
totAVeeC = 0.0;
totAVeeX = 0.0;
totAVel = 0.0;
totADI = 0.0;

totASelfT = 0.0;
totASelfVne = 0.0;
totASelfVeeC = 0.0;
totASelfVeeX = 0.0;
totASelfVel = 0.0;

# B(Complex, A)
totBVne = 0.0;
totBVen = 0.0;
totBVnn = 0.0;
totBVeeC = 0.0;
totBVeeX = 0.0;
totBVel = 0.0;
totBDI = 0.0;

totBSelfT = 0.0;
totBSelfVne = 0.0;
totBSelfVeeC = 0.0;
totBSelfVeeX = 0.0;
totBSelfVel = 0.0;

# A(Complex, Free)
totAFreeVne = 0.0;
totAFreeVen = 0.0;
totAFreeVnn = 0.0;
totAFreeVeeC = 0.0;
totAFreeVeeX = 0.0;
totAFreeVel = 0.0;
totAFreeDI = 0.0;

totAFreeSelfT = 0.0;
totAFreeSelfVne = 0.0;
totAFreeSelfVeeC = 0.0;
totAFreeSelfVeeX = 0.0;
totAFreeSelfVel = 0.0;

# B(Complex, Free)
totBFreeVne = 0.0;
totBFreeVen = 0.0;
totBFreeVnn = 0.0;
totBFreeVeeC = 0.0;
totBFreeVeeX = 0.0;
totBFreeVel = 0.0;
totBFreeDI = 0.0;

totBFreeSelfT = 0.0;
totBFreeSelfVne = 0.0;
totBFreeSelfVeeC = 0.0;
totBFreeSelfVeeX = 0.0;
totBFreeSelfVel = 0.0;

# A(Opt, Free)
totAOptVne = 0.0;
totAOptVen = 0.0;
totAOptVnn = 0.0;
totAOptVeeC = 0.0;
totAOptVeeX = 0.0;
totAOptVel = 0.0;
totAOptDI = 0.0;

totAOptSelfT = 0.0;
totAOptSelfVne = 0.0;
totAOptSelfVeeC = 0.0;
totAOptSelfVeeX = 0.0;
totAOptSelfVel = 0.0;

# B(Opt, Free)
totBOptVne = 0.0;
totBOptVen = 0.0;
totBOptVnn = 0.0;
totBOptVeeC = 0.0;
totBOptVeeX = 0.0;
totBOptVel = 0.0;
totBOptDI = 0.0;

totBOptSelfT = 0.0;
totBOptSelfVne = 0.0;
totBOptSelfVeeC = 0.0;
totBOptSelfVeeX = 0.0;
totBOptSelfVel = 0.0;

printf("\n");

# Nomenclature
printf(" Nomenclature:\n");
printf(" Fragment(Geometry,Vicinity)\n");
printf(" Geometry - Complex, or Opt\n");
printf(" Vicinity - A, B, or Free\n");
printf("\n");

# Print final output, calculate sums
# AB Interaction
printf(" A(Complex,B) Interaction with B(Complex,A)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA     Vne(atomA,fragB)     Ven(atomA,fragB)     Vnn(atomA,fragB)    VeeC(atomA,fragB)    VeeX(atomA,fragB)     Vel(atomA,fragB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NA; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", A[i], ABVne[i], ABVen[i], ABVnn[i], ABVeeC[i], ABVeeX[i], ABVel[i]);
    totABVne += ABVne[i];
    totABVen += ABVen[i];
    totABVnn += ABVnn[i];
    totABVeeC += ABVeeC[i];
    totABVeeX += ABVeeX[i];
    totABVel += ABVel[i];    
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totABVne, totABVen, totABVnn, totABVeeC, totABVeeX, totABVel);
printf("\n");
printf("\n");
printf("\n");

# A(Complex,B) Self Interaction
printf(" A(Complex,B) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA     Vne(atomA,fragA)     Ven(atomA,fragA)     Vnn(atomA,fragA)    VeeC(atomA,fragA)    VeeX(atomA,fragA)     Vel(atomA,fragA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NA; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", A[i], AVne[i], AVen[i], AVnn[i], AVeeC[i], AVeeX[i], AVel[i]);
    totAVne += AVne[i];
    totAVen += AVen[i];
    totAVnn += AVnn[i];
    totAVeeC += AVeeC[i];
    totAVeeX += AVeeX[i];
    totAVel += AVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totAVne/2, totAVen/2, totAVnn/2, totAVeeC/2, totAVeeX/2, totAVel/2);
printf("\n");

# A(Complex,B) Self Energy
printf(" A(Complex,B) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA          T(atomA)            Vne(atomA)          VeeC(atomA)          VeeX(atomA)           Vel(atomA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NA; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", A[i], ASelfT[i], ASelfVne[i], ASelfVeeC[i], ASelfVeeX[i], ASelfVel[i]);
    totASelfT += ASelfT[i];
    totASelfVne += ASelfVne[i];
    totASelfVeeC += ASelfVeeC[i];
    totASelfVeeX += ASelfVeeX[i];
    totASelfVel += ASelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totASelfT, totASelfVne, totASelfVeeC, totASelfVeeX, totASelfVel); 
printf("\n");
printf("\n");

# B(Complex, A) Self Interaction
printf(" B(Complex,A) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB     Vne(atomB,fragB)     Ven(atomB,fragB)     Vnn(atomB,fragB)    VeeC(atomB,fragB)    VeeX(atomB,fragB)     Vel(atomB,fragB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NB; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", B[i], BVne[i], BVen[i], BVnn[i], BVeeC[i], BVeeX[i], BVel[i]);
    totBVne += BVne[i];
    totBVen += BVen[i];
    totBVnn += BVnn[i];
    totBVeeC += BVeeC[i];
    totBVeeX += BVeeX[i];
    totBVel += BVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBVne/2, totBVen/2, totBVnn/2, totBVeeC/2, totBVeeX/2, totBVel/2);
printf("\n");

# B(Complex,A) Self Energy
printf(" B(Complex,A) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB          T(atomB)            Vne(atomB)          VeeC(atomB)          VeeX(atomB)           Vel(atomB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NB; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", B[i], BSelfT[i], BSelfVne[i], BSelfVeeC[i], BSelfVeeX[i], BSelfVel[i]);
    totBSelfT += BSelfT[i];
    totBSelfVne += BSelfVne[i];
    totBSelfVeeC += BSelfVeeC[i];
    totBSelfVeeX += BSelfVeeX[i];
    totBSelfVel += BSelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBSelfT, totBSelfVne, totBSelfVeeC, totBSelfVeeX, totBSelfVel);
printf("\n");
printf("\n");
printf("\n");

# A(Complex,Free) Self Interaction
printf(" A(Complex,Free) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA     Vne(atomA,fragA)     Ven(atomA,fragA)     Vnn(atomA,fragA)    VeeC(atomA,fragA)    VeeX(atomA,fragA)     Vel(atomA,fragA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NAFree; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", AFree[i], AFreeVne[i], AFreeVen[i], AFreeVnn[i], AFreeVeeC[i], AFreeVeeX[i], AFreeVel[i]);
    totAFreeVne += AFreeVne[i];
    totAFreeVen += AFreeVen[i];
    totAFreeVnn += AFreeVnn[i];
    totAFreeVeeC += AFreeVeeC[i];
    totAFreeVeeX += AFreeVeeX[i];
    totAFreeVel += AFreeVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totAFreeVne/2, totAFreeVen/2, totAFreeVnn/2, totAFreeVeeC/2, totAFreeVeeX/2, totAFreeVel/2);
printf("\n");

# A(Complex,Free) Self Energy
printf(" A(Complex,Free) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA          T(atomA)            Vne(atomA)          VeeC(atomA)          VeeX(atomA)           Vel(atomA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NAFree; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", AFree[i], AFreeSelfT[i], AFreeSelfVne[i], AFreeSelfVeeC[i], AFreeSelfVeeX[i], AFreeSelfVel[i]);
    totAFreeSelfT += AFreeSelfT[i];
    totAFreeSelfVne += AFreeSelfVne[i];
    totAFreeSelfVeeC += AFreeSelfVeeC[i];
    totAFreeSelfVeeX += AFreeSelfVeeX[i];
    totAFreeSelfVel += AFreeSelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totAFreeSelfT, totAFreeSelfVne, totAFreeSelfVeeC, totAFreeSelfVeeX, totAFreeSelfVel);
printf("\n");
printf("\n");

# B(Complex,Free) Self Interaction
printf(" B(Complex,Free) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB     Vne(atomB,fragB)     Ven(atomB,fragB)     Vnn(atomB,fragB)    VeeC(atomB,fragB)    VeeX(atomB,fragB)     Vel(atomB,fragB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NBFree; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", BFree[i], BFreeVne[i], BFreeVen[i], BFreeVnn[i], BFreeVeeC[i], BFreeVeeX[i], BFreeVel[i]);
    totBFreeVne += BFreeVne[i];
    totBFreeVen += BFreeVen[i];
    totBFreeVnn += BFreeVnn[i];
    totBFreeVeeC += BFreeVeeC[i];
    totBFreeVeeX += BFreeVeeX[i];
    totBFreeVel += BFreeVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBFreeVne/2, totBFreeVen/2, totBFreeVnn/2, totBFreeVeeC/2, totBFreeVeeX/2, totBFreeVel/2);
printf("\n");

# B(Complex,Free) Self Energy
printf(" B(Complex,Free) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB          T(atomB)              Vne(atomB)          VeeC(atomB)          VeeX(atomB)           Vel(atomB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NBFree; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", BFree[i], BFreeSelfT[i], BFreeSelfVne[i], BFreeSelfVeeC[i], BFreeSelfVeeX[i], BFreeSelfVel[i]);
    totBFreeSelfT += BFreeSelfT[i];
    totBFreeSelfVne += BFreeSelfVne[i];
    totBFreeSelfVeeC += BFreeSelfVeeC[i];
    totBFreeSelfVeeX += BFreeSelfVeeX[i];
    totBFreeSelfVel += BFreeSelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBFreeSelfT, totBFreeSelfVne, totBFreeSelfVeeC, totBFreeSelfVeeX, totBFreeSelfVel);
printf("\n");
printf("\n");
printf("\n");

# A(Opt,Free) Self Interaction
printf(" A(Opt,Free) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA     Vne(atomA,fragA)     Ven(atomA,fragA)     Vnn(atomA,fragA)    VeeC(atomA,fragA)    VeeX(atomA,fragA)     Vel(atomA,fragA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NAOpt; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", AOpt[i], AOptVne[i], AOptVen[i], AOptVnn[i], AOptVeeC[i], AOptVeeX[i], AOptVel[i]);
    totAOptVne += AOptVne[i];
    totAOptVen += AOptVen[i];
    totAOptVnn += AOptVnn[i];
    totAOptVeeC += AOptVeeC[i];
    totAOptVeeX += AOptVeeX[i];
    totAOptVel += AOptVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totAOptVne/2, totAOptVen/2, totAOptVnn/2, totAOptVeeC/2, totAOptVeeX/2, totAOptVel/2);
printf("\n");

# A(Opt,Free) Self Energy
printf(" A(Opt,Free) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomA          T(atomA)            Vne(atomA)          VeeC(atomA)          VeeX(atomA)           Vel(atomA)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NAOpt; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", AOpt[i], AOptSelfT[i], AOptSelfVne[i], AOptSelfVeeC[i], AOptSelfVeeX[i], AOptSelfVel[i]);
    totAOptSelfT += AOptSelfT[i];
    totAOptSelfVne += AOptSelfVne[i];
    totAOptSelfVeeC += AOptSelfVeeC[i];
    totAOptSelfVeeX += AOptSelfVeeX[i];
    totAOptSelfVel += AOptSelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totAOptSelfT, totAOptSelfVne, totAOptSelfVeeC, totAOptSelfVeeX, totAOptSelfVel);
printf("\n");
printf("\n");

# B(Opt, Free) Self Interaction
printf(" B(Opt,Free) Self Interaction\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB     Vne(atomB,fragB)     Ven(atomB,fragB)     Vnn(atomB,fragB)    VeeC(atomB,fragB)    VeeX(atomB,fragB)     Vel(atomB,fragB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NBOpt; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", BOpt[i], BOptVne[i], BOptVen[i], BOptVnn[i], BOptVeeC[i], BOptVeeX[i], BOptVel[i]);
    totBOptVne += BOptVne[i];
    totBOptVen += BOptVen[i];
    totBOptVnn += BOptVnn[i];
    totBOptVeeC += BOptVeeC[i];
    totBOptVeeX += BOptVeeX[i];
    totBOptVel += BOptVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum/2%20.10f %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBOptVne/2, totBOptVen/2, totBOptVnn/2, totBOptVeeC/2, totBOptVeeX/2, totBOptVel/2);
printf("\n");

# B(Opt,Free) Self Energy
printf(" B(Opt,Free) Self Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" atomB          T(atomB)              Vne(atomB)          VeeC(atomB)          VeeX(atomB)           Vel(atomB)\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
for(i=1; i<NBOpt; i++) {
    printf("%5s %20.10f %20.10f %20.10f %20.10f %20.10f\n", BOpt[i], BOptSelfT[i], BOptSelfVne[i], BOptSelfVeeC[i], BOptSelfVeeX[i], BOptSelfVel[i]);
    totBOptSelfT += BOptSelfT[i];
    totBOptSelfVne += BOptSelfVne[i];
    totBOptSelfVeeC += BOptSelfVeeC[i];
    totBOptSelfVeeX += BOptSelfVeeX[i];
    totBOptSelfVel += BOptSelfVel[i];
}
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum  %20.10f %20.10f %20.10f %20.10f %20.10f\n", totBOptSelfT, totBOptSelfVne, totBOptSelfVeeC, totBOptSelfVeeX, totBOptSelfVel);
printf("\n");
printf("\n");

# Calculate components of electronic deformation energy
# A
AElDefT = totASelfT - totAFreeSelfT;

AIntElDefVeeX = totAVeeX/2 - totAFreeVeeX/2;
ASelfElDefVeeX = totASelfVeeX - totAFreeSelfVeeX;
AElDefVeeX = AIntElDefVeeX + ASelfElDefVeeX;

AIntElDefVel = totAVel/2 - totAFreeVel/2;
ASelfElDefVel = totASelfVel - totAFreeSelfVel;
AElDefVel = AIntElDefVel + ASelfElDefVel;

# B
BElDefT = totBSelfT - totBFreeSelfT;

BIntElDefVeeX = totBVeeX/2 - totBFreeVeeX/2;
BSelfElDefVeeX = totBSelfVeeX - totBFreeSelfVeeX;
BElDefVeeX = BIntElDefVeeX + BSelfElDefVeeX;

BIntElDefVel = totBVel/2 - totBFreeVel/2;
BSelfElDefVel = totBSelfVel - totBFreeSelfVel;
BElDefVel = BIntElDefVel + BSelfElDefVel;

# Calculate components of geometry deformation energy
# A
AGeoDefT = totAFreeSelfT - totAOptSelfT;

AIntGeoDefVeeX = totAFreeVeeX/2 - totAOptVeeX/2;
ASelfGeoDefVeeX = totAFreeSelfVeeX - totAOptSelfVeeX;
AGeoDefVeeX = AIntGeoDefVeeX + ASelfGeoDefVeeX;

AIntGeoDefVel = totAFreeVel/2 - totAOptVel/2;
ASelfGeoDefVel = totAFreeSelfVel - totAOptSelfVel;
AGeoDefVel = AIntGeoDefVel + ASelfGeoDefVel;

# B
BGeoDefT = totBFreeSelfT - totBOptSelfT;

BIntGeoDefVeeX = totBFreeVeeX/2 - totBOptVeeX/2;
BSelfGeoDefVeeX = totBFreeSelfVeeX - totBOptSelfVeeX;
BGeoDefVeeX = BIntGeoDefVeeX + BSelfGeoDefVeeX;

BIntGeoDefVel = totBFreeVel/2 - totBOptVel/2;
BSelfGeoDefVel = totBFreeSelfVel - totBOptSelfVel;
BGeoDefVel = BIntGeoDefVel + BSelfGeoDefVel;

# Tot
AGeoDef = AGeoDefT + AGeoDefVeeX + AGeoDefVel;
BGeoDef = BGeoDefT + BGeoDefVeeX + BGeoDefVel;
ABGeoDef = AGeoDef + BGeoDef;

# Calculate final sums
sumAB = totABVeeX + totABVel;

# El
sumAIntElDef = AIntElDefVel + AIntElDefVeeX;
sumASelfElDef = AElDefT + ASelfElDefVel + ASelfElDefVeeX;
sumAElDef = sumAIntElDef + sumASelfElDef;

sumBIntElDef = BIntElDefVel + BIntElDefVeeX;
sumBSelfElDef = BElDefT + BSelfElDefVel + BSelfElDefVeeX;
sumBElDef = sumBIntElDef + sumBSelfElDef;

sumABElDef = sumAElDef + sumBElDef;

sumElDefT = AElDefT + BElDefT;
sumElDefVeeX = AElDefVeeX + BElDefVeeX;
sumElDefVel = AElDefVel + BElDefVel;

# Geo
sumAIntGeoDef = AIntGeoDefVel + AIntGeoDefVeeX;
sumASelfGeoDef = AGeoDefT + ASelfGeoDefVel + ASelfGeoDefVeeX;
sumAGeoDef = sumAIntGeoDef + sumASelfGeoDef;

sumBIntGeoDef = BIntGeoDefVel + BIntGeoDefVeeX;
sumBSelfGeoDef = BGeoDefT + BSelfGeoDefVel + BSelfGeoDefVeeX;
sumBGeoDef = sumBIntGeoDef + sumBSelfGeoDef;

sumABGeoDef = sumAGeoDef + sumBGeoDef;

sumGeoDefT = AGeoDefT + BGeoDefT;
sumGeoDefVeeX = AGeoDefVeeX + BGeoDefVeeX;
sumGeoDefVel = AGeoDefVel + BGeoDefVel;

sumDefT = sumElDefT + sumGeoDefT;
sumDefVeeX = totABVeeX + sumElDefVeeX + sumGeoDefVeeX;
sumDefVel = totABVel + sumElDefVel + sumGeoDefVel;

sumFinal = sumAB + sumABElDef + sumABGeoDef;

# Interaction energy
printf("\n");
printf(" AB Total Interaction Energy\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf("      Component                                             VeeX [kcal/mol]       Vel [kcal/mol]      Sum\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf("         AB                             %30.2f %20.2f %13.2f\n", totABVeeX*kcal, totABVel*kcal, sumAB*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");

# Electronic Deformation energy Free -> Complex
printf("\n");
printf(" A and B Electronic Deformation Energies\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf("             Component                    T [kcal/mol]      VeeX [kcal/mol]       Vel [kcal/mol]     Sum\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Intra: A(Complex,Free) -> A(Complex,B)%10.2f %20.2f %20.2f %13.2f\n", AElDefT*kcal, ASelfElDefVeeX*kcal, ASelfElDefVel*kcal, sumASelfElDef*kcal);
printf(" Inter: A(Complex,Free) -> A(Complex,B)%10s %20.2f %20.2f %13.2f\n", "-----", AIntElDefVeeX*kcal, AIntElDefVel*kcal, sumAIntElDef*kcal);
printf(" Tot:   A(Complex,Free) -> A(Complex,B)%10.2f %20.2f %20.2f %13.2f\n", AElDefT*kcal, AElDefVeeX*kcal, AElDefVel*kcal, sumAElDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Intra: B(Complex,Free) -> B(Complex,A)%10.2f %20.2f %20.2f %13.2f\n", BElDefT*kcal, BSelfElDefVeeX*kcal, BSelfElDefVel*kcal, sumBSelfElDef*kcal);
printf(" Inter: B(Complex,Free) -> B(Complex,A)%10s %20.2f %20.2f %13.2f\n", "-----", BIntElDefVeeX*kcal, BIntElDefVel*kcal, sumBIntElDef*kcal);
printf(" Tot:   B(Complex,Free) -> B(Complex,A)%10.2f %20.2f %20.2f %13.2f\n", BElDefT*kcal, BElDefVeeX*kcal, BElDefVel*kcal, sumBElDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Final                       %20.2f %20.2f %20.2f %13.2f\n", sumElDefT*kcal, sumElDefVeeX*kcal, sumElDefVel*kcal, sumABElDef*kcal);

# Geometry Deformation energy Optimized -> Free
printf("\n");
printf(" A and B Geometry Deformation Energies\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf("             Component                    T [kcal/mol]      VeeX [kcal/mol]       Vel [kcal/mol]     Sum\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Intra: A(Opt,Free) -> A(Complex,Free) %10.2f %20.2f %20.2f %13.2f\n", AGeoDefT*kcal, ASelfGeoDefVeeX*kcal, ASelfGeoDefVel*kcal, sumASelfGeoDef*kcal);
printf(" Inter: A(Opt,Free) -> A(Complex,Free) %10s %20.2f %20.2f %13.2f\n", "-----", AIntGeoDefVeeX*kcal, AIntGeoDefVel*kcal, sumAIntGeoDef*kcal);
printf(" Tot:   A(Opt,Free) -> A(Complex,Free) %10.2f %20.2f %20.2f %13.2f\n", AGeoDefT*kcal, AGeoDefVeeX*kcal, AGeoDefVel*kcal, sumAGeoDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Intra: B(Opt,Free) -> B(Complex,Free) %10.2f %20.2f %20.2f %13.2f\n", BGeoDefT*kcal, BSelfGeoDefVeeX*kcal, BSelfGeoDefVel*kcal, sumBSelfGeoDef*kcal);
printf(" Inter: B(Opt,Free) -> B(Complex,Free) %10s %20.2f %20.2f %13.2f\n", "-----", BIntGeoDefVeeX*kcal, BIntGeoDefVel*kcal, sumBIntGeoDef*kcal);
printf(" Tot:   B(Opt,Free) -> B(Complex,Free) %10.2f %20.2f %20.2f %13.2f\n", BGeoDefT*kcal, BGeoDefVeeX*kcal, BGeoDefVel*kcal, sumBGeoDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Final                       %20.2f %20.2f %20.2f %13.2f\n", sumGeoDefT*kcal, sumGeoDefVeeX*kcal, sumGeoDefVel*kcal, sumABGeoDef*kcal);

# Final components of binding energy
printf("\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf("      Component                           T [kcal/mol]      VeeX [kcal/mol]       Vel [kcal/mol]       Sum\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Interaction energy     %46.2f %20.2f %14.2f\n",                           totABVeeX*kcal, totABVel*kcal, sumAB*kcal);
printf(" Electronic deformation %25.2f %20.2f %20.2f %14.2f\n", sumElDefT*kcal, sumElDefVeeX*kcal, sumElDefVel*kcal, sumABElDef*kcal);
printf(" Geometry deformation   %25.2f %20.2f %20.2f %14.2f\n", sumGeoDefT*kcal, sumGeoDefVeeX*kcal, sumGeoDefVel*kcal, sumABGeoDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" Sum                    %25.2f %20.2f %20.2f %14.2f\n", sumDefT*kcal, sumDefVeeX*kcal, sumDefVel*kcal, sumFinal*kcal );

# Final results
printf("\n");
printf("\n");
printf(" FINAL RESULTS\n");
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" INTERACTION ENERGY              %16.2f kcal/mol\n", sumAB*kcal);
printf(" ELECTRONIC DEFORMATION ENERGY   %16.2f kcal/mol\n", sumABElDef*kcal);
printf(" GEOMETRY DEFORMATION ENERGY     %16.2f kcal/mol\n", sumABGeoDef*kcal);
printf("-------------------------------------------------------------------------------------------------------------------------------------\n");
printf(" FINAL IQA BINDING ENERGY        %16.2f kcal/mol\n", sumFinal*kcal);

# Print xyz files for control of correct selection of FragmentA and FragmentB
if( control == 1 ) {
    printf("\n");
    printf("Control output is switched ON\n");
    printf("Writing:\n")

    # Complex AB
    printf("== AB.xyz\n")
    print N-1 > "AB.xyz";
    print "" >> "AB.xyz";
    for(i=1; i<N; i++) {
        element = AB[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, ABx[i]*au, ABy[i]*au, ABz[i]*au) >> "AB.xyz";
    }

    # Fragment A
    printf("== A.xyz\n")
    print NA-1 > "A.xyz";
    print "" >> "A.xyz";
    for(i=1; i<NA; i++) {
        element = A[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Ax[i]*au, Ay[i]*au, Az[i]*au) >> "A.xyz";
    }

    # Fragment B
    printf("== B.xyz\n")
    print NB-1 > "B.xyz";
    print "" >> "B.xyz";
    for(i=1; i<NB; i++) {
        element = B[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, Bx[i]*au, By[i]*au, Bz[i]*au) >> "B.xyz";
    }

    # Reference A - Free
    printf("== AFree.xyz\n")
    print NAFree-1 > "AFree.xyz";
    print "" >> "AFree.xyz";
    for(i = 1; i < NAFree; i++) {
        element = AFree[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, AFreex[i]*au, AFreey[i]*au, AFreez[i]*au) >> "AFree.xyz";
    }

    # Reference B - Free
    printf("== BFree.xyz\n")
    print NBFree-1 > "BFree.xyz";
    print "" >> "BFree.xyz";
    for(i = 1; i < NBFree; i++) {
        element = BFree[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, BFreex[i]*au, BFreey[i]*au, BFreez[i]*au) >> "BFree.xyz";
    }

    # Reference A - Opt
    printf("== AOpt.xyz\n")
    print NAOpt-1 > "AOpt.xyz";
    print "" >> "AOpt.xyz";
    for(i = 1; i < NAOpt; i++) {
        element = AOpt[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, AOptx[i]*au, AOpty[i]*au, AOptz[i]*au) >> "AOpt.xyz";
    }

    # Reference B - Opt
    printf("== BOpt.xyz\n")
    print NBOpt-1 > "BOpt.xyz";
    print "" >> "BOpt.xyz";
    for(i = 1; i < NBOpt; i++) {
        element = BOpt[i];
        sub(/[0-9]{1,}/, "", element);
        printf("%s\t%f\t%f\t%f\n", element, BOptx[i]*au, BOpty[i]*au, BOptz[i]*au) >> "BOpt.xyz";
    }
}

}
