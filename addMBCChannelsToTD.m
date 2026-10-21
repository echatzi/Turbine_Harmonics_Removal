function TD = addMBCChannelsToTD(TD)
% Add the six Coleman MBC channels to an already-loaded TD structure.

if ~isstruct(TD) || ~isscalar(TD)
    error('TD must be a scalar structure containing channel fields.');
end

gearRatio = 126.793;
bladeOrder = {'A', 'C', 'B'};
[tDay, genRpm] = getSignal(TD, 'GenSdSpd');
tSec = (tDay - tDay(1)) * 86400;
if numel(tDay) < 2 || any(~isfinite(tSec)) || any(diff(tSec) <= 0)
    error('GenSdSpd must contain at least two strictly increasing timestamps.');
end
if any(~isfinite(genRpm))
    error('GenSdSpd contains non-finite values; azimuth cannot be integrated.');
end

bladeFlp = zeros(numel(tSec), 3);
bladeEdg = zeros(numel(tSec), 3);
for bladeIndex = 1:3
    bladeId = bladeOrder{bladeIndex};
    xLE = loadAligned(TD, sprintf('BLSSMV_Bld%s_LE', bladeId), tDay, tSec);
    xTE = loadAligned(TD, sprintf('BLSSMV_Bld%s_TE', bladeId), tDay, tSec);
    xPS = loadAligned(TD, sprintf('BLSSMV_Bld%s_PS', bladeId), tDay, tSec);
    xSS = loadAligned(TD, sprintf('BLSSMV_Bld%s_SS', bladeId), tDay, tSec);
    bladeFlp(:, bladeIndex) = xPS - xSS;
    bladeEdg(:, bladeIndex) = xLE - xTE;
end

azimuth = cumtrapz(tSec, max(genRpm(:) / gearRatio, 0) * (2*pi/60));
if isfield(TD, 'Orientat_Azimuth') && isfield(TD.Orientat_Azimuth, 'value')
    measuredAzimuth = loadAligned(TD, 'Orientat_Azimuth', tDay, tSec);
    firstValidAzimuth = find(isfinite(measuredAzimuth), 1);
    if ~isempty(firstValidAzimuth)
        azimuth = azimuth + deg2rad(measuredAzimuth(firstValidAzimuth));
    end
end

thetaA = azimuth;
thetaC = azimuth + 2*pi/3;
thetaB = azimuth + 4*pi/3;
[Flp_Q0, Flp_Qc, Flp_Qs] = mbcForward( ...
    bladeFlp(:,1), bladeFlp(:,2), bladeFlp(:,3), thetaA, thetaC, thetaB);
[Edg_Q0, Edg_Qc, Edg_Qs] = mbcForward( ...
    bladeEdg(:,1), bladeEdg(:,2), bladeEdg(:,3), thetaA, thetaC, thetaB);

channelNames = {'Flp_Q0', 'Flp_Qc', 'Flp_Qs', 'Edg_Q0', 'Edg_Qc', 'Edg_Qs'};
channelValues = {Flp_Q0, Flp_Qc, Flp_Qs, Edg_Q0, Edg_Qc, Edg_Qs};
for channelIndex = 1:numel(channelNames)
    TD.(channelNames{channelIndex}) = struct( ...
        'time', tDay(:), 'value', channelValues{channelIndex}(:));
end
end

function [t, x] = getSignal(TD, name)
if ~isfield(TD, name) || ~isfield(TD.(name), 'time') || ~isfield(TD.(name), 'value')
    error('Required channel is missing or malformed: %s', name);
end
t = TD.(name).time(:);
x = TD.(name).value(:);
if numel(t) ~= numel(x)
    error('Channel time/value lengths differ: %s', name);
end
end

function y = loadAligned(TD, channelName, tDayRef, tSecRef)
if ~isfield(TD, channelName)
    error('Required channel is missing: %s', channelName);
end
[tChannel, xChannel] = getSignal(TD, channelName);
tChannelSec = (tChannel - tDayRef(1)) * 86400;
y = interp1(tChannelSec, xChannel, tSecRef, 'linear', NaN);
end

function [Q0, Qc, Qs] = mbcForward(QA, QB, QC, thetaA, thetaB, thetaC)
Q0 = (QA + QB + QC) / 3;
Qc = (2/3) * (QA.*cos(thetaA) + QB.*cos(thetaB) + QC.*cos(thetaC));
Qs = (2/3) * (QA.*sin(thetaA) + QB.*sin(thetaB) + QC.*sin(thetaC));
end