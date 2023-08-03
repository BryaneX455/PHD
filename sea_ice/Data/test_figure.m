floe_idx = [ 1,   38,   49,   52, 1610, 1611, 1621, 1622, ...
          1625, 1626, 1627, 1628, 1630, 1631, 1635, 1639, ...
          1642, 1643, 1654, 1656, 1657, 1659, 1660, 1661, ...
          1663, 1669, 1670, 1673, 1677, 1692, 1712, 1760, ...
          1837, 1960, 1963, 1975, 1977, 1980, 2000];

gprdata = load("GprInterpData.mat");
lindata = load("LinearInterpData.mat");

figure
hold on
for iIdx = 1:length(floe_idx)
    plot(gprdata.xInterpArray(floe_idx(iIdx), :), gprdata.yInterpArray(floe_idx(iIdx), :), 'r')
    plot(lindata.xInterpArray(floe_idx(iIdx), :), lindata.yInterpArray(floe_idx(iIdx), :), 'k--')
end