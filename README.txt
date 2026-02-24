This dataset consists of data captured by MATLAB's HDL Verifier tool. It represents ultrasonic amplitudes, converted to voltages by transducing receivers, and quantized by the ADC on the Zybo Z7 FPGA. In this folder one may also find the code required to process these data.

Each data folder is broken into multiple sub-folders. Because there were 16 transmission-mode receivers, multiplexing the receivers was necessary. ctrl_00, _01, _10, and _11 represent the different control bits sent to the multiplexer, and therefore a unique set of 4 receivers. The rx subfolder is for the 4 reflection receivers.

Within each data subfolder is a .mat file containing either reflection or transmission data for a particular beam angle. The .mat file itself is a struct containing receiver amplitudes over time, time markers for the beginnings of transmitted bursts, and the beam angle. This struct is unpacked in the get_receiver_data function.

Using the code is straightforward. Under %%% FOLDER PATH INPUTS %%%, one can specify the folder path for reflection mode receivers (folderPath) and the folder path for transmission mode receivers (folderPathTx), along with their respective calibration data folder paths. A few key results can be replicated with the following paths and parameters:

function G = dspglobals
    G.INTERPOLATION_TYPE = 'nearest';
    G.VOLTAGE_AUTORANGE = true;
    G.MAX_VOLTAGE = 475;
    G.MIN_VOLTAGE = 0;
    G.TRANSMISSION_VOLTAGE_AUTORANGE = true;
    G.TRANSMISSION_RECEIVERS_MAX_VOLTAGE = 360;
    G.TRANSMISSION_RECEIVERS_MIN_VOLTAGE = 0;
    G.PROMINENCE_THRESHOLD_MULTIPLIER = 0.3;
    G.VOLTAGE_THRESHOLD_MULTIPLIER = 0.2;
    G.DO_X_MASK = false;
    G.YMASK_TOLERANCE = 3;
    G.XMASK_TOLERANCE = 30;
end

1. Plastic and cotton on the right:

folderPath = './11-12 plastic_cotton right 3 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-12 plastic_cotton right 3 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';

2. Plastic and cotton in the center:

folderPath = './11-8 plastic_cotton center 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-8 plastic_cotton center 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';

3. Cotton and plastic on the left:

folderPath = './11-8 plastic_cotton left 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-8 plastic_cotton left 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';

4. Cotton alone on the right:

folderPath = './11-12 cotton right 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-12 cotton right 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';
G.VOLTAGE_AUTORANGE = false;
G.TRANSMISSION_VOLTAGE_AUTORANGE = false;

5. Cotton alone on the left:

folderPath = './11-8 cotton left 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-8 cotton left 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';
G.VOLTAGE_AUTORANGE = false;
G.TRANSMISSION_VOLTAGE_AUTORANGE = false;

6. Cotton alone in the middle:

folderPath = './11-8 cotton center 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './11-8 cotton center 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';
G.VOLTAGE_AUTORANGE = false;
G.TRANSMISSION_VOLTAGE_AUTORANGE = false;

(Decreasing the value of PROMINENCE_THRESHOLD_MULTIPLIER can help highlight the distinct peaks for these next two.)

7. Two plastics, separated laterally by 11mm:

folderPath = './2-11 1p1_36p3 cm_sep 15sweep 30v/rx/';
folderPathCali = './2-8 air 15sweep 30v/rx/';
folderPathTx = './2-11 1p1_36p3 cm_sep 15sweep 30v/';
folderPathTxCali = './2-8 air 15sweep 30v/';
G.INTERPOLATION_TYPE = 'natural';

8. Two plastics, separated axially by 103 mm:

folderPath = './3-1 10p3_31p1 2 cm_axsep 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-1 10p3_31p1 2 cm_axsep 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';
G.INTERPOLATION_TYPE = 'natural';

9. Two plastics, separated axially by 121 mm:

folderPath = './3-1 12p1_31p1 cm_axsep 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-1 12p1_31p1 cm_axsep 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

10. Two plastics separated:

folderPath = './1-8 2plastics 5cm 3 15sweep 30v/rx/';
folderPathCali = './11-23 air 15sweep 30v/rx/';
folderPathTx = './1-8 2plastics 5cm 3 15sweep 30v/';
folderPathTxCali = './11-23 air 15sweep 30v/';

(Table 1 data):

11. 5x5 cm plastic:

folderPath = './3-25 cotton_plastic_5mm 4 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_plastic_5mm 4 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

12. 4x4 cm plastic:

folderPath = './3-25 cotton_plastic_4mm 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_plastic_4mm 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

13. 3x3 cm plastic:

folderPath = './3-25 cotton_plastic_3mm 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_plastic_3mm 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

14. 2x2 cm plastic:

folderPath = './3-25 cotton_plastic_2mm 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_plastic_2mm 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

15. 1x1 cm plastic:

folderPath = './3-25 cotton_plastic_1mm 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_plastic_1mm 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

16. Cotton:

folderPath = './3-25 cotton 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';

17. Seeded cotton:

folderPath = './3-25 cotton_seeded 15sweep 30v/rx/';
folderPathCali = './3-11 air 15sweep 30v/rx/';
folderPathTx = './3-25 cotton_seeded 15sweep 30v/';
folderPathTxCali = './3-11 air 15sweep 30v/';


More data examples are available upon request.