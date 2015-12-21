{
Anomaly=0.75*exp(- (pow(X-250,2) + pow(Y-160,2) + pow(Z-190,2)) /800);
Model=Background+Anomaly;
if(Topocheck==0 ) {Model=-100;

}
}

