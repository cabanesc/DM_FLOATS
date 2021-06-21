function myebcc(X,Y,ERR,coly,colerr);
%
% myeb(Y,varargin);
%
% This function makes nice coloured, shaded error bars. Exactly what
% it does depends on Y, and on whether you give it one or two inputs. 
%
% If you only pass it Y, and no other arguments, it assuemd you're
% giving it raw data. 
%
%		myeb(Raw_Data)
%
% 	.) if Y is 2D array, it will then plot mean(Y) with errorbars given
% 	by std(Y). In this case there is only one mean vector with its
% 	errorbars. 
% 
%	.) if Y is 3D array, it will plot size(Y,3) lines with the
%	associated errorbars. Line k will be mean(Y(:,:,k)) with errorbars
%	given by std(Y(:,:,k))
%
% If you pass it 2 arguments, each has to be at most 2D. 
%
%		myeb(mu,std)
%
% 	.) if mu and std are 1D, it just plots one line given by mu with a
% 	shaded region given by std. 
%
%	.) if mu and std are 2D, it will plot size(Y,2) lines in the
%	standard sequence of colours; each line mu(:,k) will have a shaded
%	region in the same colour, but less saturated given by std(:,k)
%
%
% Quentin Huys, 2007
% Center for Theoretical Neuroscience, Columbia University
% Email: qhuys [at] n e u r o theory [dot] columbia.edu
% (just get rid of the spaces, replace [at] with @ and [dot] with .)



colerr_edge=colerr-0.2;

	m=Y;
	s=ERR;

        ind1=[1:length(X)];
	ind2=fliplr(ind1);
	
	hold on; 
	h=fill([X fliplr(X)],[m-s m(ind2)+s(ind2)],colerr);
	set(h,'edgecolor',colerr_edge)
	plot(X,m,'linewidth',2,'Color',coly)
	