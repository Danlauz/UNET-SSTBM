dir='C:\Users\dany-\OneDrive\Bureau\Article7_UNet_PumpingTest\128x128_BC_S_Q\Image300T60V';

times=[0 1 5 15 30 60 120 180 720]*60;

load('DataTrainValidation.mat')

nbsim=size(Data,4);
nbtrain=nbsim*0.8;
nbValid=nbsim*0.16;
nbTest=nbsim*0.04;
nx=size(Data,1); ny=size(Data,2);

b=3;
n=0.3; %Aquifer porosity
cW=4.4*10^-10; % Water compressibility
rhow=1000; % Water density
g=9.81; %gravity
Ss=rhow*g*(cR+n*cW);
S=Ss*b;

%% Building images for the CNN training

a=0;b=0;
for j=1:size(Data,3)-2
    jj=1;
    for i=1:nbsim
        a=a+1;
        % Create floating point image.
        im=Data(:,:,j+2,i);
        rgbImage = cat(3,im);
        % Image must be single precision.
        rgbImage = single(rgbImage);
        % Create tiff object.
        if i<=nbtrain
            fileName = [dir '\OutputTiff\Training\' sprintf('R%d.tif',a)];
        elseif (nbtrain<i) && (i<=nbtrain+nbValid)
            fileName =[dir '\OutputTiff\Validation\' sprintf('R%d.tif',a)];
        elseif (nbtrain+nbValid)<i
            fileName = [dir '\OutputTiff\Test\' sprintf('R%d.tif',a)];
        end
        t = Tiff(fileName, 'w');
        % Set tags.
        tagstruct.ImageLength = size(rgbImage, 1);
        tagstruct.ImageWidth = size(rgbImage, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.LinearRaw;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = 1;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        t.setTag(tagstruct);
        t.write(double(rgbImage));
        t.close();
        % Recall image.
        m2 = imread(fileName);
        % Check that it's the same as what we wrote out.
        maxDiff = max(max(m2-rgbImage));
    end

    for i=1:nbsim
        b=b+1;
        % Create floating point image.
        im1=  Data(:,:,1,i);
        im2=  Data(:,:,2,i);
        im3=  Data(:,:,j+1,i);      
        im4=zeros(nx,ny); im4(nx/2,ny/2)=times(j);
        if j==1
            im5=zeros(nx,ny);
        else
            im5=zeros(nx,ny); im5(nx/2,ny/2)=Q(i)*minute/litre;
        end
        im6=zeros(nx,ny); im6(nx/2,ny/2)=log10(Ss(i));
        rgbImage = cat(3,im1,im2,im3,im4,im5,im6);
        % Image must be single precision.
        rgbImage = single(rgbImage);
        % Create tiff object.
        if i<=nbtrain
            fileName = [dir '\InputTiff\Training\' sprintf('R%d.tif',b)];
        elseif (nbtrain<i) && (i<=nbtrain+nbValid)
            fileName = [dir '\InputTiff\Validation\' sprintf('R%d.tif',b)];
        elseif (nbtrain+nbValid)<i
            fileName = [dir '\InputTiff\Test\' sprintf('R%d.tif',b)];
        end
        t = Tiff(fileName, 'w');
        % Set tags.
        tagstruct.ImageLength = size(rgbImage, 1);
        tagstruct.ImageWidth = size(rgbImage, 2);
        tagstruct.Compression = Tiff.Compression.None;
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.Photometric = Tiff.Photometric.LinearRaw;
        tagstruct.BitsPerSample = 64;
        tagstruct.SamplesPerPixel = 6;
        tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        t.setTag(tagstruct);
        t.write(double(rgbImage));
        t.close();
        % Recall image.
        m3 = imread(fileName);
        % Check that it's the same as what we wrote out.
        maxDiff = max(max(m3-rgbImage));
    end
end