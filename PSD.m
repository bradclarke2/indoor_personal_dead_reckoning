function [px]=PSD(x,showplot)
%[px]=jonsPSD(x,showplot) returns unweighted power spectrum for data with time values in column 1 and magnitude values in column 2
% If 'showplot' is 'y' a plot of the power spectrum will be produced

if nargin < 2
    showplot='n';
end

tinc=x(2,1)-x(1,1);
Gx=fft(x(:,2));

%Compute power spectrum
pxx=(((1/length(Gx))^2)*Gx.*conj(Gx));                            %the squared magnitude of the transform is obtained from the conjugate of the complex fx  
                                                                %and is scaled by dividing by the number of points (also squared)
px(1,2)=pxx(1);                                                 %the first data point is the dc level and is placed in first column of px
                                                                %Transform is a mirror image about half the sample freq, it needs to be folded in half
px(2:0.5*length(Gx),2)=2*(pxx(2:0.5*length(Gx)));               %so multiply by 2 and divide by sqrt(2) to get rms level. Because we're working with power spectrum, needs to be (2/sqrt(2))^2

px(:,1)=1/(length(Gx)*tinc)*((1:(0.5*length(Gx)))-1)';          %frequency values obtained knowing frequency separation between points = 1/(no. of samples * time increment)

if showplot=='y'
%     figure
%     semilogy(px(:,1),px(:,2))
% %     xlim ([0 5000])
    figure
    plot(px(:,1),px(:,2))
 
end


