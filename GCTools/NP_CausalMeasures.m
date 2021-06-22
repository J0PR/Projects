function [DTF,PDC,DC,GGC,GCGeweke,PartitionsGCGeweke,freqStruct] = NP_CausalMeasures(X,Fs,NLAGS)

    [H Z S np_freqs]=calcHandCovNParam(X,Fs,0);
    
    %np_GCtest=freqGCtest(X,'np',H,Z,S,np_freqs,Fs);
    multiWaitbar( 'Calculating measures.', 0.2);
    GCGeweke=GewekeGC1(X,'np',H,Z,S,np_freqs,Fs);
    multiWaitbar( 'Calculating measures.', 0.4);
    PartitionsGCGeweke=GewekeGC1Partitions(X,'np',H,Z,S,np_freqs,Fs);
    multiWaitbar( 'Calculating measures.', 0.6);
    %np_PartitionsGCPCA=PCAGCPartitions(X,'np',H,Z,S,np_freqs,Fs);
    multiWaitbar( 'Calculating measures.', 0.8);
    [S,h,PDC,DTF,DC,GGC]=mvfreqz('np',H,Z,np_freqs,Fs);
    multiWaitbar( 'Calculating measures.', 1);
    multiWaitbar( 'Calculating measures.', 'Close');
    freqStruct=np_freqs;
    
%    Written by João Rodrigues
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.