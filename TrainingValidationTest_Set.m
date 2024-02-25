%% Training and Validation sets for CNN

dataFolder1 = fullfile('C:\Users\dany-\OneDrive\Bureau\Article7_UNet_PumpingTest\128x128_BC_S_Q\Image500T100V\InputTiff\Training');
dataFolder2 = fullfile('C:\Users\dany-\OneDrive\Bureau\Article7_UNet_PumpingTest\128x128_BC_S_Q\Image500T100V\OutputTiff\Training');

imds1 = imageDatastore(dataFolder1,"IncludeSubfolders",true,"FileExtensions",".tif");
imds2 = imageDatastore(dataFolder2,"IncludeSubfolders",true,"FileExtensions",".tif");

imdsTrain_CNN=combine(imds1,imds2);


dataFolder1 = fullfile('C:\Users\dany-\OneDrive\Bureau\Article7_UNet_PumpingTest\128x128_BC_S_Q\Image500T100V\InputTiff\Validation');
dataFolder2 = fullfile('C:\Users\dany-\OneDrive\Bureau\Article7_UNet_PumpingTest\128x128_BC_S_Q\Image500T100V\OutputTiff\Validation');

imds1 = imageDatastore(dataFolder1,"IncludeSubfolders",true,"FileExtensions",".tif");
imds2 = imageDatastore(dataFolder2,"IncludeSubfolders",true,"FileExtensions",".tif");

imdsValidation_CNN=combine(imds1,imds2);


