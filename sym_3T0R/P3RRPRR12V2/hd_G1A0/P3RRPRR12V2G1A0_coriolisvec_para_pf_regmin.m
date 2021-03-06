% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRPRR12V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:15:54
% EndTime: 2020-08-06 19:16:19
% DurationCPUTime: 25.24s
% Computational Cost: add. (191368->738), mult. (262065->1371), div. (11064->12), fcn. (135963->18), ass. (0->562)
t1478 = xDP(2);
t1479 = xDP(1);
t1480 = pkin(5) - pkin(6);
t1412 = pkin(1) * t1478 + t1479 * t1480;
t1468 = sin(qJ(2,3));
t1816 = qJ(3,3) * t1478;
t1355 = t1412 * t1468 + t1816;
t1413 = pkin(1) * t1479 - t1478 * t1480;
t1815 = qJ(3,3) * t1479;
t1358 = t1413 * t1468 + t1815;
t1655 = t1468 * t1815;
t1364 = t1413 + 0.2e1 * t1655;
t1656 = t1468 * t1816;
t1370 = 0.2e1 * t1656 + t1412;
t1841 = pkin(1) * t1468;
t1427 = qJ(3,3) + t1841;
t1842 = legFrame(3,3);
t1440 = sin(t1842);
t1474 = cos(qJ(2,3));
t1458 = t1474 ^ 2;
t1465 = pkin(1) * qJ(3,3);
t1469 = sin(qJ(1,3));
t1477 = xDP(3);
t1481 = pkin(2) + pkin(3);
t1845 = cos(qJ(1,3));
t1392 = t1478 * t1469 + t1479 * t1845;
t1395 = -t1469 * t1479 + t1478 * t1845;
t1669 = cos(t1842);
t1530 = t1392 * t1669 + t1395 * t1440;
t1737 = (qJ(3,3) + t1481) * (-qJ(3,3) + t1481);
t1615 = t1458 * t1737;
t1716 = t1477 * t1481;
t1681 = -0.2e1 * t1716;
t1301 = t1530 * t1615 + (-t1477 * (-t1468 * t1737 + t1465) + ((t1364 * t1845 + t1370 * t1469) * t1669 + (-t1364 * t1469 + t1370 * t1845) * t1440) * t1481) * t1474 + t1427 * t1716 + (t1681 * t1458 + (t1355 * t1469 + t1358 * t1845) * t1669 + (t1355 * t1845 - t1358 * t1469) * t1440) * qJ(3,3);
t1818 = qJ(3,3) * t1468;
t1422 = pkin(1) + t1818;
t1722 = t1474 * t1481;
t1406 = t1422 + t1722;
t1397 = 0.1e1 / t1406;
t1483 = 0.1e1 / qJ(3,3);
t1759 = t1397 * t1483;
t1298 = t1301 * t1759;
t1365 = t1413 + t1655;
t1371 = t1656 + t1412;
t1322 = (-qJ(3,3) * t1477 + t1481 * t1530) * t1458 + ((t1365 * t1845 + t1371 * t1469) * t1669 + (-t1365 * t1469 + t1371 * t1845) * t1440 + t1468 * t1716) * t1474 + t1477 * t1427;
t1645 = t1322 * t1759;
t1316 = pkin(2) * t1645;
t1274 = t1316 - t1298;
t1265 = -pkin(3) * t1645 - t1274;
t1353 = -t1440 * t1392 + t1395 * t1669;
t1728 = t1468 * t1480;
t1623 = t1353 * t1728;
t1714 = t1481 * t1483;
t1314 = (-t1322 * t1714 + t1623) * t1397;
t1457 = t1474 * t1458;
t1482 = qJ(3,3) ^ 2;
t1769 = t1353 * t1397;
t1569 = t1737 * t1769;
t1657 = qJ(3,3) * t1769;
t1598 = t1458 * t1657;
t1727 = t1468 * t1483;
t1642 = t1322 * t1727;
t1398 = 0.1e1 / t1406 ^ 2;
t1484 = 0.1e1 / qJ(3,3) ^ 2;
t1757 = t1398 * t1484;
t1643 = t1322 * t1757;
t1688 = 0.2e1 * t1481;
t1493 = pkin(1) ^ 2;
t1690 = t1480 ^ 2 + t1493;
t1758 = t1398 * t1483;
t1868 = t1480 * t1481;
t1295 = t1480 * t1298;
t1340 = pkin(1) * t1657;
t1881 = -t1295 - 0.2e1 * t1340;
t1236 = -(-t1457 * t1569 + (-t1353 * t1422 * t1688 - t1322 * t1480) * t1397 * t1458 + (t1881 * t1468 + (-t1353 * (t1482 + t1690) + t1642 * t1868) * t1397) * t1474) * t1353 * t1758 - (-t1480 * t1598 + ((-pkin(3) * t1322 * t1483 + t1623) * t1397 - t1274) * t1722 + t1265 * t1422) * t1643 - (-t1314 * t1474 + t1422 * t1645) * t1301 * t1757;
t1837 = pkin(2) * t1236;
t1794 = t1265 * t1481;
t1253 = -qJ(3,3) * t1322 * t1397 + t1794;
t1689 = pkin(5) ^ 2 + t1493;
t1417 = (-0.2e1 * pkin(5) + pkin(6)) * pkin(6) + t1689;
t1437 = 0.2e1 * t1458 - 0.1e1;
t1557 = 0.2e1 * t1868;
t1625 = t1481 * t1769;
t1761 = t1397 * t1468;
t1626 = t1353 * t1761;
t1766 = t1353 * t1483;
t1767 = t1353 * t1480;
t1852 = -0.3e1 * t1482;
t1859 = t1481 ^ 2;
t1867 = t1482 - t1859;
t1878 = 0.2e1 * pkin(1);
t1849 = -pkin(5) / 0.2e1;
t1885 = -0.2e1 * t1481 * (t1849 + pkin(6) / 0.2e1);
t1882 = (t1484 * (((t1253 * t1481 + t1569 * t1728) * t1474 + pkin(1) * t1794 + (t1253 * t1468 + (-t1353 * t1437 * t1868 - pkin(1) * t1322) * t1397) * qJ(3,3)) * t1322 + (-t1314 * t1722 + ((pkin(1) * t1483 + t1468) * t1481 * t1322 + (t1458 - 0.1e1) * qJ(3,3) * t1767) * t1397) * t1301) + (-(t1859 + t1852) * t1457 * t1625 + (-0.3e1 * (-t1482 / 0.3e1 + t1859) * t1818 + t1867 * t1878) * t1769 * t1458 + (((-t1295 - 0.4e1 * t1340) * t1481 - t1867 * t1480 * t1645) * t1468 + (t1852 - t1417) * t1625) * t1474 + ((-t1557 * t1645 + t1295) * t1458 + (-t1417 - t1482) * t1626 + t1645 * t1885 + t1881) * qJ(3,3)) * t1766) * t1398;
t1888 = t1882 + t1837;
t1470 = sin(qJ(2,2));
t1820 = qJ(3,2) * t1478;
t1356 = t1412 * t1470 + t1820;
t1819 = qJ(3,2) * t1479;
t1359 = t1413 * t1470 + t1819;
t1658 = t1470 * t1819;
t1366 = t1413 + 0.2e1 * t1658;
t1659 = t1470 * t1820;
t1372 = 0.2e1 * t1659 + t1412;
t1840 = pkin(1) * t1470;
t1428 = qJ(3,2) + t1840;
t1843 = legFrame(2,3);
t1441 = sin(t1843);
t1475 = cos(qJ(2,2));
t1460 = t1475 ^ 2;
t1466 = pkin(1) * qJ(3,2);
t1471 = sin(qJ(1,2));
t1846 = cos(qJ(1,2));
t1391 = -t1471 * t1479 + t1846 * t1478;
t1393 = t1478 * t1471 + t1479 * t1846;
t1670 = cos(t1843);
t1531 = t1391 * t1441 + t1393 * t1670;
t1736 = (qJ(3,2) + t1481) * (-qJ(3,2) + t1481);
t1614 = t1460 * t1736;
t1302 = t1531 * t1614 + (-t1477 * (-t1470 * t1736 + t1466) + ((t1366 * t1846 + t1372 * t1471) * t1670 + (-t1366 * t1471 + t1372 * t1846) * t1441) * t1481) * t1475 + t1428 * t1716 + (t1681 * t1460 + (t1356 * t1471 + t1359 * t1846) * t1670 + (t1356 * t1846 - t1359 * t1471) * t1441) * qJ(3,2);
t1822 = qJ(3,2) * t1470;
t1424 = pkin(1) + t1822;
t1720 = t1475 * t1481;
t1407 = t1424 + t1720;
t1400 = 0.1e1 / t1407;
t1486 = 0.1e1 / qJ(3,2);
t1754 = t1400 * t1486;
t1299 = t1302 * t1754;
t1367 = t1413 + t1658;
t1373 = t1659 + t1412;
t1323 = (-qJ(3,2) * t1477 + t1481 * t1531) * t1460 + ((t1367 * t1846 + t1373 * t1471) * t1670 + (-t1367 * t1471 + t1373 * t1846) * t1441 + t1470 * t1716) * t1475 + t1477 * t1428;
t1641 = t1323 * t1754;
t1317 = pkin(2) * t1641;
t1276 = t1317 - t1299;
t1266 = -pkin(3) * t1641 - t1276;
t1352 = t1391 * t1670 - t1441 * t1393;
t1726 = t1470 * t1480;
t1627 = t1352 * t1726;
t1713 = t1481 * t1486;
t1313 = (-t1323 * t1713 + t1627) * t1400;
t1459 = t1475 * t1460;
t1485 = qJ(3,2) ^ 2;
t1773 = t1352 * t1400;
t1570 = t1736 * t1773;
t1660 = qJ(3,2) * t1773;
t1599 = t1460 * t1660;
t1725 = t1470 * t1486;
t1638 = t1323 * t1725;
t1401 = 0.1e1 / t1407 ^ 2;
t1487 = 0.1e1 / qJ(3,2) ^ 2;
t1752 = t1401 * t1487;
t1639 = t1323 * t1752;
t1753 = t1401 * t1486;
t1296 = t1480 * t1299;
t1341 = pkin(1) * t1660;
t1880 = -t1296 - 0.2e1 * t1341;
t1235 = -(-t1459 * t1570 + (-t1352 * t1424 * t1688 - t1323 * t1480) * t1400 * t1460 + (t1880 * t1470 + (-t1352 * (t1485 + t1690) + t1638 * t1868) * t1400) * t1475) * t1352 * t1753 - (-t1480 * t1599 + ((-pkin(3) * t1323 * t1486 + t1627) * t1400 - t1276) * t1720 + t1266 * t1424) * t1639 - (-t1313 * t1475 + t1424 * t1641) * t1302 * t1752;
t1838 = pkin(2) * t1235;
t1793 = t1266 * t1481;
t1254 = -qJ(3,2) * t1323 * t1400 + t1793;
t1438 = 0.2e1 * t1460 - 0.1e1;
t1629 = t1481 * t1773;
t1756 = t1400 * t1470;
t1630 = t1352 * t1756;
t1770 = t1352 * t1486;
t1771 = t1352 * t1480;
t1851 = -0.3e1 * t1485;
t1866 = t1485 - t1859;
t1883 = (t1487 * (((t1254 * t1481 + t1570 * t1726) * t1475 + pkin(1) * t1793 + (t1254 * t1470 + (-t1352 * t1438 * t1868 - pkin(1) * t1323) * t1400) * qJ(3,2)) * t1323 + (-t1313 * t1720 + ((pkin(1) * t1486 + t1470) * t1481 * t1323 + (t1460 - 0.1e1) * qJ(3,2) * t1771) * t1400) * t1302) + (-(t1859 + t1851) * t1459 * t1629 + (-0.3e1 * (-t1485 / 0.3e1 + t1859) * t1822 + t1866 * t1878) * t1773 * t1460 + (((-t1296 - 0.4e1 * t1341) * t1481 - t1866 * t1480 * t1641) * t1470 + (t1851 - t1417) * t1629) * t1475 + ((-t1557 * t1641 + t1296) * t1460 + (-t1417 - t1485) * t1630 + t1641 * t1885 + t1880) * qJ(3,2)) * t1770) * t1401;
t1887 = t1883 + t1838;
t1472 = sin(qJ(2,1));
t1824 = qJ(3,1) * t1478;
t1357 = t1412 * t1472 + t1824;
t1823 = qJ(3,1) * t1479;
t1360 = t1413 * t1472 + t1823;
t1661 = t1472 * t1823;
t1368 = t1413 + 0.2e1 * t1661;
t1662 = t1472 * t1824;
t1374 = 0.2e1 * t1662 + t1412;
t1839 = pkin(1) * t1472;
t1429 = qJ(3,1) + t1839;
t1844 = legFrame(1,3);
t1442 = sin(t1844);
t1476 = cos(qJ(2,1));
t1462 = t1476 ^ 2;
t1467 = pkin(1) * qJ(3,1);
t1473 = sin(qJ(1,1));
t1847 = cos(qJ(1,1));
t1394 = t1478 * t1473 + t1847 * t1479;
t1396 = -t1473 * t1479 + t1478 * t1847;
t1671 = cos(t1844);
t1532 = t1394 * t1671 + t1396 * t1442;
t1735 = (qJ(3,1) + t1481) * (-qJ(3,1) + t1481);
t1613 = t1462 * t1735;
t1303 = t1532 * t1613 + (-t1477 * (-t1472 * t1735 + t1467) + ((t1368 * t1847 + t1374 * t1473) * t1671 + (-t1368 * t1473 + t1374 * t1847) * t1442) * t1481) * t1476 + t1429 * t1716 + (t1681 * t1462 + (t1357 * t1473 + t1360 * t1847) * t1671 + (t1357 * t1847 - t1360 * t1473) * t1442) * qJ(3,1);
t1826 = qJ(3,1) * t1472;
t1426 = pkin(1) + t1826;
t1718 = t1476 * t1481;
t1408 = t1426 + t1718;
t1403 = 0.1e1 / t1408;
t1489 = 0.1e1 / qJ(3,1);
t1749 = t1403 * t1489;
t1300 = t1303 * t1749;
t1369 = t1413 + t1661;
t1375 = t1662 + t1412;
t1324 = (-qJ(3,1) * t1477 + t1481 * t1532) * t1462 + ((t1369 * t1847 + t1375 * t1473) * t1671 + (-t1369 * t1473 + t1375 * t1847) * t1442 + t1472 * t1716) * t1476 + t1477 * t1429;
t1637 = t1324 * t1749;
t1318 = pkin(2) * t1637;
t1278 = t1318 - t1300;
t1267 = -pkin(3) * t1637 - t1278;
t1354 = -t1442 * t1394 + t1396 * t1671;
t1724 = t1472 * t1480;
t1619 = t1354 * t1724;
t1712 = t1481 * t1489;
t1315 = (-t1324 * t1712 + t1619) * t1403;
t1461 = t1476 * t1462;
t1488 = qJ(3,1) ^ 2;
t1765 = t1354 * t1403;
t1568 = t1735 * t1765;
t1663 = qJ(3,1) * t1765;
t1600 = t1462 * t1663;
t1723 = t1472 * t1489;
t1634 = t1324 * t1723;
t1404 = 0.1e1 / t1408 ^ 2;
t1490 = 0.1e1 / qJ(3,1) ^ 2;
t1747 = t1404 * t1490;
t1635 = t1324 * t1747;
t1748 = t1404 * t1489;
t1297 = t1480 * t1300;
t1342 = pkin(1) * t1663;
t1879 = -t1297 - 0.2e1 * t1342;
t1237 = -(-t1461 * t1568 + (-t1354 * t1426 * t1688 - t1324 * t1480) * t1403 * t1462 + (t1879 * t1472 + (-t1354 * (t1488 + t1690) + t1634 * t1868) * t1403) * t1476) * t1354 * t1748 - (-t1480 * t1600 + ((-pkin(3) * t1324 * t1489 + t1619) * t1403 - t1278) * t1718 + t1267 * t1426) * t1635 - (-t1315 * t1476 + t1426 * t1637) * t1303 * t1747;
t1836 = pkin(2) * t1237;
t1792 = t1267 * t1481;
t1255 = -qJ(3,1) * t1324 * t1403 + t1792;
t1439 = 0.2e1 * t1462 - 0.1e1;
t1621 = t1481 * t1765;
t1751 = t1403 * t1472;
t1622 = t1354 * t1751;
t1762 = t1354 * t1489;
t1763 = t1354 * t1480;
t1850 = -0.3e1 * t1488;
t1865 = t1488 - t1859;
t1884 = (t1490 * (((t1255 * t1481 + t1568 * t1724) * t1476 + pkin(1) * t1792 + (t1255 * t1472 + (-t1354 * t1439 * t1868 - pkin(1) * t1324) * t1403) * qJ(3,1)) * t1324 + (-t1315 * t1718 + ((pkin(1) * t1489 + t1472) * t1481 * t1324 + (t1462 - 0.1e1) * qJ(3,1) * t1763) * t1403) * t1303) + (-(t1859 + t1850) * t1461 * t1621 + (-0.3e1 * (-t1488 / 0.3e1 + t1859) * t1826 + t1865 * t1878) * t1765 * t1462 + (((-t1297 - 0.4e1 * t1342) * t1481 - t1865 * t1480 * t1637) * t1472 + (t1850 - t1417) * t1621) * t1476 + ((-t1557 * t1637 + t1297) * t1462 + (-t1417 - t1488) * t1622 + t1637 * t1885 + t1879) * qJ(3,1)) * t1762) * t1404;
t1886 = t1884 + t1836;
t1271 = t1316 - t1298 / 0.2e1;
t1877 = -0.4e1 * t1271;
t1272 = t1317 - t1299 / 0.2e1;
t1876 = -0.4e1 * t1272;
t1273 = t1318 - t1300 / 0.2e1;
t1875 = -0.4e1 * t1273;
t1644 = t1322 * t1758;
t1606 = 0.2e1 * t1644;
t1874 = t1397 * t1606;
t1640 = t1323 * t1753;
t1605 = 0.2e1 * t1640;
t1873 = t1400 * t1605;
t1636 = t1324 * t1748;
t1604 = 0.2e1 * t1636;
t1872 = t1403 * t1604;
t1832 = pkin(2) * t1474;
t1871 = t1468 * (pkin(1) + t1832);
t1831 = pkin(2) * t1475;
t1870 = t1470 * (pkin(1) + t1831);
t1830 = pkin(2) * t1476;
t1869 = t1472 * (pkin(1) + t1830);
t1861 = -0.2e1 * pkin(1);
t1858 = -0.2e1 * t1274;
t1857 = -0.2e1 * t1276;
t1856 = -0.2e1 * t1278;
t1855 = 0.4e1 * t1458;
t1854 = 0.4e1 * t1460;
t1853 = 0.4e1 * t1462;
t1848 = pkin(5) / 0.2e1;
t1835 = pkin(2) * t1458;
t1834 = pkin(2) * t1460;
t1833 = pkin(2) * t1462;
t1817 = qJ(3,3) * t1474;
t1409 = t1468 * t1481 - t1817;
t1399 = t1397 * t1398;
t1624 = t1399 * t1766;
t1588 = t1322 * t1624;
t1768 = t1353 * t1398;
t1247 = -(t1298 * t1468 + (t1767 + (-t1468 * t1714 + t1474) * t1322) * t1397) * t1768 + t1409 * t1588 - t1468 * t1301 * t1624;
t1829 = pkin(5) * t1247;
t1821 = qJ(3,2) * t1475;
t1410 = t1470 * t1481 - t1821;
t1402 = t1400 * t1401;
t1628 = t1402 * t1770;
t1587 = t1323 * t1628;
t1772 = t1352 * t1401;
t1248 = -(t1299 * t1470 + (t1771 + (-t1470 * t1713 + t1475) * t1323) * t1400) * t1772 + t1410 * t1587 - t1470 * t1302 * t1628;
t1828 = pkin(5) * t1248;
t1825 = qJ(3,1) * t1476;
t1411 = t1472 * t1481 - t1825;
t1405 = t1403 * t1404;
t1620 = t1405 * t1762;
t1586 = t1324 * t1620;
t1764 = t1354 * t1404;
t1249 = -(t1300 * t1472 + (t1763 + (-t1472 * t1712 + t1476) * t1324) * t1403) * t1764 + t1411 * t1586 - t1472 * t1303 * t1620;
t1827 = pkin(5) * t1249;
t1814 = t1458 * qJ(3,3);
t1813 = t1460 * qJ(3,2);
t1812 = t1462 * qJ(3,1);
t1810 = 0.2e1 * pkin(2);
t1809 = t1235 * t1475;
t1808 = t1236 * t1474;
t1807 = t1237 * t1476;
t1806 = t1247 * t1397;
t1805 = t1247 * t1468;
t1804 = t1247 * t1474;
t1803 = t1247 * t1483;
t1802 = t1248 * t1400;
t1801 = t1248 * t1470;
t1800 = t1248 * t1475;
t1799 = t1248 * t1486;
t1798 = t1249 * t1403;
t1797 = t1249 * t1472;
t1796 = t1249 * t1476;
t1795 = t1249 * t1489;
t1319 = t1322 ^ 2;
t1791 = t1319 * t1484;
t1320 = t1323 ^ 2;
t1790 = t1320 * t1487;
t1321 = t1324 ^ 2;
t1789 = t1321 * t1490;
t1433 = t1480 * t1845;
t1361 = t1469 * t1427 - t1433 * t1468;
t1376 = t1440 * t1469 - t1669 * t1845;
t1685 = 0.2e1 * t1818;
t1421 = pkin(1) + t1685;
t1382 = t1421 * t1469 - t1433;
t1430 = t1469 * t1480;
t1541 = t1427 * t1845 + t1468 * t1430;
t1564 = t1421 * t1845 + t1430;
t1325 = -t1376 * t1615 - (t1440 * t1382 - t1564 * t1669) * t1722 - (t1440 * t1361 - t1541 * t1669) * qJ(3,3);
t1788 = t1325 * t1468;
t1377 = t1440 * t1845 + t1669 * t1469;
t1326 = t1377 * t1615 + (t1382 * t1669 + t1440 * t1564) * t1722 + (t1361 * t1669 + t1440 * t1541) * qJ(3,3);
t1787 = t1326 * t1468;
t1434 = t1480 * t1846;
t1362 = t1471 * t1428 - t1434 * t1470;
t1378 = t1441 * t1471 - t1670 * t1846;
t1686 = 0.2e1 * t1822;
t1423 = pkin(1) + t1686;
t1384 = t1423 * t1471 - t1434;
t1431 = t1471 * t1480;
t1540 = t1428 * t1846 + t1470 * t1431;
t1562 = t1423 * t1846 + t1431;
t1327 = -t1378 * t1614 - (t1441 * t1384 - t1562 * t1670) * t1720 - (t1441 * t1362 - t1540 * t1670) * qJ(3,2);
t1786 = t1327 * t1470;
t1379 = t1441 * t1846 + t1670 * t1471;
t1328 = t1379 * t1614 + (t1384 * t1670 + t1441 * t1562) * t1720 + (t1362 * t1670 + t1441 * t1540) * qJ(3,2);
t1785 = t1328 * t1470;
t1435 = t1480 * t1847;
t1363 = t1473 * t1429 - t1435 * t1472;
t1380 = t1442 * t1473 - t1671 * t1847;
t1687 = 0.2e1 * t1826;
t1425 = pkin(1) + t1687;
t1386 = t1425 * t1473 - t1435;
t1432 = t1473 * t1480;
t1539 = t1429 * t1847 + t1472 * t1432;
t1560 = t1425 * t1847 + t1432;
t1329 = -t1380 * t1613 - (t1442 * t1386 - t1560 * t1671) * t1718 - (t1442 * t1363 - t1539 * t1671) * qJ(3,1);
t1784 = t1329 * t1472;
t1381 = t1442 * t1847 + t1671 * t1473;
t1330 = t1381 * t1613 + (t1386 * t1671 + t1442 * t1560) * t1718 + (t1363 * t1671 + t1442 * t1539) * qJ(3,1);
t1783 = t1330 * t1472;
t1383 = t1422 * t1469 - t1433;
t1563 = t1422 * t1845 + t1430;
t1334 = t1376 * t1722 + t1383 * t1440 - t1563 * t1669;
t1782 = t1334 * t1474;
t1335 = t1377 * t1722 + t1383 * t1669 + t1440 * t1563;
t1781 = t1335 * t1474;
t1385 = t1424 * t1471 - t1434;
t1561 = t1424 * t1846 + t1431;
t1336 = t1378 * t1720 + t1385 * t1441 - t1561 * t1670;
t1780 = t1336 * t1475;
t1337 = t1379 * t1720 + t1385 * t1670 + t1441 * t1561;
t1779 = t1337 * t1475;
t1387 = t1426 * t1473 - t1435;
t1559 = t1426 * t1847 + t1432;
t1338 = t1380 * t1718 + t1387 * t1442 - t1559 * t1671;
t1778 = t1338 * t1476;
t1339 = t1381 * t1718 + t1387 * t1671 + t1442 * t1559;
t1777 = t1339 * t1476;
t1349 = t1352 ^ 2;
t1346 = t1349 * t1401;
t1776 = t1349 * t1460;
t1350 = t1353 ^ 2;
t1347 = t1350 * t1398;
t1775 = t1350 * t1458;
t1351 = t1354 ^ 2;
t1348 = t1351 * t1404;
t1774 = t1351 * t1462;
t1760 = t1397 * t1474;
t1755 = t1400 * t1475;
t1750 = t1403 * t1476;
t1414 = -pkin(2) * t1468 + t1817;
t1746 = t1414 * t1474;
t1415 = -pkin(2) * t1470 + t1821;
t1745 = t1415 * t1475;
t1416 = -pkin(2) * t1472 + t1825;
t1744 = t1416 * t1476;
t1743 = t1421 * t1474;
t1742 = t1422 * t1474;
t1741 = t1423 * t1475;
t1740 = t1424 * t1475;
t1739 = t1425 * t1476;
t1738 = t1426 * t1476;
t1454 = t1468 ^ 2;
t1734 = t1454 * t1483;
t1455 = t1470 ^ 2;
t1733 = t1455 * t1486;
t1456 = t1472 ^ 2;
t1732 = t1456 * t1489;
t1731 = t1458 * t1483;
t1730 = t1460 * t1486;
t1729 = t1462 * t1489;
t1721 = t1474 * t1483;
t1719 = t1475 * t1486;
t1717 = t1476 * t1489;
t1275 = t1316 - 0.2e1 * t1298;
t1492 = pkin(2) ^ 2;
t1443 = -t1482 + t1492;
t1542 = qJ(3,3) * t1808 - t1888 * t1468;
t1609 = t1769 * t1858;
t1618 = t1397 * t1721;
t1711 = -t1542 * pkin(5) - (pkin(1) * t1685 + t1443 * t1458 + t1742 * t1810 + t1482 + t1689) * t1247 + t1598 * t1877 - 0.2e1 * ((t1275 * t1849 + t1340) * t1322 - (-pkin(2) * t1301 + t1322 * t1443) * t1626) * t1618 - (-pkin(5) * t1319 * t1758 + pkin(1) * t1609) * t1468 - qJ(3,3) * t1609;
t1277 = t1317 - 0.2e1 * t1299;
t1444 = -t1485 + t1492;
t1543 = qJ(3,2) * t1809 - t1887 * t1470;
t1608 = t1773 * t1857;
t1617 = t1400 * t1719;
t1710 = -t1543 * pkin(5) - (pkin(1) * t1686 + t1444 * t1460 + t1740 * t1810 + t1485 + t1689) * t1248 + t1599 * t1876 - 0.2e1 * ((t1277 * t1849 + t1341) * t1323 - (-pkin(2) * t1302 + t1323 * t1444) * t1630) * t1617 - (-pkin(5) * t1320 * t1753 + pkin(1) * t1608) * t1470 - qJ(3,2) * t1608;
t1279 = t1318 - 0.2e1 * t1300;
t1445 = -t1488 + t1492;
t1544 = qJ(3,1) * t1807 - t1886 * t1472;
t1607 = t1765 * t1856;
t1616 = t1403 * t1717;
t1709 = -t1544 * pkin(5) - (pkin(1) * t1687 + t1445 * t1462 + t1738 * t1810 + t1488 + t1689) * t1249 + t1600 * t1875 - 0.2e1 * (t1324 * (t1279 * t1849 + t1342) - (-pkin(2) * t1303 + t1324 * t1445) * t1622) * t1616 - (-pkin(5) * t1321 * t1748 + pkin(1) * t1607) * t1472 - qJ(3,1) * t1607;
t1684 = -0.2e1 * t1814;
t1692 = pkin(2) * t1347 + 0.2e1 * t1301 * t1643;
t1708 = (t1492 + t1482) * t1236 + t1414 * t1829 + pkin(2) * t1882 + t1692 * qJ(3,3) + (-(-t1443 * t1468 + t1465) * t1474 + (t1684 + t1841) * pkin(2)) * t1347;
t1683 = -0.2e1 * t1813;
t1693 = pkin(2) * t1346 + 0.2e1 * t1302 * t1639;
t1707 = (t1492 + t1485) * t1235 + t1415 * t1828 + pkin(2) * t1883 + t1693 * qJ(3,2) + (-(-t1444 * t1470 + t1466) * t1475 + (t1683 + t1840) * pkin(2)) * t1346;
t1682 = -0.2e1 * t1812;
t1691 = pkin(2) * t1348 + 0.2e1 * t1303 * t1635;
t1706 = (t1492 + t1488) * t1237 + t1416 * t1827 + pkin(2) * t1884 + t1691 * qJ(3,1) + (-(-t1445 * t1472 + t1467) * t1476 + (t1682 + t1839) * pkin(2)) * t1348;
t1651 = t1319 * t1757;
t1674 = pkin(5) * t1805;
t1705 = t1674 - qJ(3,3) * (t1347 + t1651) + (t1814 - t1871) * t1347 - t1888;
t1649 = t1320 * t1752;
t1673 = pkin(5) * t1801;
t1704 = t1673 - qJ(3,2) * (t1346 + t1649) + (t1813 - t1870) * t1346 - t1887;
t1647 = t1321 * t1747;
t1672 = pkin(5) * t1797;
t1703 = t1672 - qJ(3,1) * (t1348 + t1647) + (t1812 - t1869) * t1348 - t1886;
t1603 = pkin(5) * t1651;
t1677 = pkin(5) * t1808;
t1702 = -t1677 - (0.2e1 * t1871 + (-0.2e1 * t1458 + 0.2e1) * qJ(3,3)) * t1247 + t1468 * t1603 - (t1606 * t1743 + (t1271 * t1855 + t1858) * t1397) * t1353;
t1602 = pkin(5) * t1649;
t1679 = pkin(5) * t1809;
t1701 = -t1679 - (0.2e1 * t1870 + (-0.2e1 * t1460 + 0.2e1) * qJ(3,2)) * t1248 + t1470 * t1602 - (t1605 * t1741 + (t1272 * t1854 + t1857) * t1400) * t1352;
t1601 = pkin(5) * t1647;
t1675 = pkin(5) * t1807;
t1700 = -t1675 - (0.2e1 * t1869 + (-0.2e1 * t1462 + 0.2e1) * qJ(3,1)) * t1249 + t1472 * t1601 - (t1604 * t1739 + (t1273 * t1853 + t1856) * t1403) * t1354;
t1678 = pkin(5) * t1236 * t1468;
t1699 = t1678 - 0.2e1 * (t1742 + t1835) * t1247 - (t1626 * t1877 - t1603) * t1474 - (-0.2e1 * t1427 * t1483 + t1855) * t1322 * t1768;
t1680 = pkin(5) * t1235 * t1470;
t1698 = t1680 - 0.2e1 * (t1740 + t1834) * t1248 - (t1630 * t1876 - t1602) * t1475 - (-0.2e1 * t1428 * t1486 + t1854) * t1323 * t1772;
t1676 = pkin(5) * t1237 * t1472;
t1697 = t1676 - 0.2e1 * (t1738 + t1833) * t1249 - (t1622 * t1875 - t1601) * t1476 - (-0.2e1 * t1429 * t1489 + t1853) * t1324 * t1764;
t1696 = pkin(5) * t1800 + 0.2e1 * qJ(3,2) * t1235 + (-t1741 - 0.2e1 * t1834) * t1346 + t1693;
t1695 = pkin(5) * t1804 + 0.2e1 * qJ(3,3) * t1236 + (-t1743 - 0.2e1 * t1835) * t1347 + t1692;
t1694 = pkin(5) * t1796 + 0.2e1 * qJ(3,1) * t1237 + (-t1739 - 0.2e1 * t1833) * t1348 + t1691;
t1650 = t1399 * t1791;
t1648 = t1402 * t1790;
t1646 = t1405 * t1789;
t1633 = t1349 * t1719;
t1632 = t1350 * t1721;
t1631 = t1351 * t1717;
t1612 = t1468 * t1721;
t1611 = t1470 * t1719;
t1610 = t1472 * t1717;
t1597 = t1731 * t1806;
t1596 = t1247 * t1612;
t1595 = t1730 * t1802;
t1594 = t1248 * t1611;
t1593 = t1729 * t1798;
t1592 = t1249 * t1610;
t1591 = (t1275 * t1474 + t1322 * t1761) * t1644;
t1590 = (t1277 * t1475 + t1323 * t1756) * t1640;
t1589 = (t1279 * t1476 + t1324 * t1751) * t1636;
t1585 = t1334 * t1618;
t1584 = t1335 * t1618;
t1583 = t1336 * t1617;
t1582 = t1337 * t1617;
t1581 = t1338 * t1616;
t1580 = t1339 * t1616;
t1579 = t1401 * t1633;
t1578 = t1349 * t1402 * t1730;
t1577 = t1402 * t1633;
t1576 = t1398 * t1632;
t1575 = t1350 * t1399 * t1731;
t1574 = t1399 * t1632;
t1573 = t1404 * t1631;
t1572 = t1351 * t1405 * t1729;
t1571 = t1405 * t1631;
t1567 = t1399 * t1612;
t1566 = t1402 * t1611;
t1565 = t1405 * t1610;
t1556 = t1397 * t1596;
t1555 = t1400 * t1594;
t1554 = t1403 * t1592;
t1553 = t1470 * t1578;
t1552 = t1349 * t1566;
t1551 = t1468 * t1575;
t1550 = t1350 * t1567;
t1549 = t1472 * t1572;
t1548 = t1351 * t1565;
t1538 = -(t1804 * t1878 - t1678) * t1397 + (t1322 * t1721 * t1848 + t1353 * t1841) * t1874;
t1537 = -(t1800 * t1878 - t1680) * t1400 + (t1323 * t1719 * t1848 + t1352 * t1840) * t1873;
t1536 = -(t1796 * t1878 - t1676) * t1403 + (t1324 * t1717 * t1848 + t1354 * t1839) * t1872;
t1535 = -(t1805 * t1861 - t1677) * t1397 + (t1474 * pkin(1) * t1353 + t1642 * t1849) * t1874;
t1534 = -(t1801 * t1861 - t1679) * t1400 + (t1475 * pkin(1) * t1352 + t1638 * t1849) * t1873;
t1533 = -(t1797 * t1861 - t1675) * t1403 + (t1476 * pkin(1) * t1354 + t1634 * t1849) * t1872;
t1529 = -t1235 * t1756 - t1475 * t1648;
t1528 = -t1235 * t1755 + t1470 * t1648;
t1527 = -t1236 * t1761 - t1474 * t1650;
t1526 = -t1236 * t1760 + t1468 * t1650;
t1525 = -t1237 * t1751 - t1476 * t1646;
t1524 = -t1237 * t1750 + t1472 * t1646;
t1523 = -0.2e1 * t1322 * t1353 * t1567 - t1454 * t1806;
t1522 = -0.2e1 * t1323 * t1352 * t1566 - t1455 * t1802;
t1521 = -0.2e1 * t1324 * t1354 * t1565 - t1456 * t1798;
t1520 = t1247 * t1734 + t1248 * t1733 + t1249 * t1732;
t1519 = -0.2e1 * t1437 * t1588 - 0.2e1 * t1760 * t1805;
t1518 = -0.2e1 * t1438 * t1587 - 0.2e1 * t1755 * t1801;
t1517 = -0.2e1 * t1439 * t1586 - 0.2e1 * t1750 * t1797;
t1510 = t1592 + t1594 + t1596;
t1509 = t1334 * t1597 + t1336 * t1595 + t1338 * t1593;
t1508 = t1335 * t1597 + t1337 * t1595 + t1339 * t1593;
t1507 = t1334 * t1556 + t1336 * t1555 + t1338 * t1554;
t1506 = t1335 * t1556 + t1337 * t1555 + t1339 * t1554;
t1390 = (pkin(1) + 0.2e1 * t1830) * t1472 + t1682 + qJ(3,1);
t1389 = (pkin(1) + 0.2e1 * t1831) * t1470 + t1683 + qJ(3,2);
t1388 = (pkin(1) + 0.2e1 * t1832) * t1468 + t1684 + qJ(3,3);
t1333 = -0.2e1 * t1404 * t1774 + t1348;
t1332 = -0.2e1 * t1401 * t1776 + t1346;
t1331 = -0.2e1 * t1398 * t1775 + t1347;
t1306 = -t1348 + (t1774 - t1789) * t1404;
t1305 = -t1346 + (t1776 - t1790) * t1401;
t1304 = -t1347 + (t1775 - t1791) * t1398;
t1213 = t1884 - t1672 + 0.2e1 * t1836;
t1212 = t1882 - t1674 + 0.2e1 * t1837;
t1211 = t1883 - t1673 + 0.2e1 * t1838;
t1207 = t1544 + 0.2e1 * t1827;
t1206 = t1542 + 0.2e1 * t1829;
t1205 = t1543 + 0.2e1 * t1828;
t1 = [-t1377 * t1806 - t1379 * t1802 - t1381 * t1798, 0, 0, t1334 * t1551 + t1336 * t1553 + t1338 * t1549 + t1377 * t1523 + t1379 * t1522 + t1381 * t1521, -t1331 * t1585 - t1332 * t1583 - t1333 * t1581 + t1377 * t1519 + t1379 * t1518 + t1381 * t1517, t1377 * t1527 + t1379 * t1529 + t1381 * t1525 - t1507, t1377 * t1526 + t1379 * t1528 + t1381 * t1524 - t1509, -t1235 * t1583 - t1236 * t1585 - t1237 * t1581, t1536 * t1381 + t1537 * t1379 + t1538 * t1377 + t1507 * pkin(5) + (-t1334 * t1550 - t1336 * t1552 - t1338 * t1548) * pkin(1), t1533 * t1381 + t1534 * t1379 + t1535 * t1377 + t1509 * pkin(5) + (-t1334 * t1575 - t1336 * t1578 - t1338 * t1572) * pkin(1), (-t1338 * t1390 - t1784) * t1571 + (-t1334 * t1388 - t1788) * t1574 + (-t1336 * t1389 - t1786) * t1577 + ((-t1213 * t1778 - t1237 * t1329) * t1489 + t1697 * t1381) * t1403 + ((-t1211 * t1780 - t1235 * t1327) * t1486 + t1698 * t1379) * t1400 + ((-t1212 * t1782 - t1236 * t1325) * t1483 + t1699 * t1377) * t1397, t1377 * t1591 + t1379 * t1590 + t1381 * t1589 + (-t1381 * t1207 + (-t1338 * t1744 + t1784) * t1795) * t1403 + (-t1379 * t1205 + (-t1336 * t1745 + t1786) * t1799) * t1400 + (-t1377 * t1206 + (-t1334 * t1746 + t1788) * t1803) * t1397, (t1700 * t1381 + (t1306 * t1329 - t1694 * t1778) * t1489) * t1403 + (t1701 * t1379 + (t1305 * t1327 - t1696 * t1780) * t1486) * t1400 + (t1702 * t1377 + (t1304 * t1325 - t1695 * t1782) * t1483) * t1397, (t1709 * t1381 + (t1329 * t1703 - t1706 * t1778) * t1489) * t1403 + (t1710 * t1379 + (t1327 * t1704 - t1707 * t1780) * t1486) * t1400 + (t1711 * t1377 + (t1325 * t1705 - t1708 * t1782) * t1483) * t1397, 0; -t1376 * t1806 - t1378 * t1802 - t1380 * t1798, 0, 0, -t1335 * t1551 - t1337 * t1553 - t1339 * t1549 + t1376 * t1523 + t1378 * t1522 + t1380 * t1521, t1331 * t1584 + t1332 * t1582 + t1333 * t1580 + t1376 * t1519 + t1378 * t1518 + t1380 * t1517, t1376 * t1527 + t1378 * t1529 + t1380 * t1525 + t1506, t1376 * t1526 + t1378 * t1528 + t1380 * t1524 + t1508, t1235 * t1582 + t1236 * t1584 + t1237 * t1580, t1536 * t1380 + t1537 * t1378 + t1538 * t1376 - t1506 * pkin(5) + (t1335 * t1550 + t1337 * t1552 + t1339 * t1548) * pkin(1), t1533 * t1380 + t1534 * t1378 + t1535 * t1376 - t1508 * pkin(5) + (t1335 * t1575 + t1337 * t1578 + t1339 * t1572) * pkin(1), (t1339 * t1390 - t1783) * t1571 + (t1335 * t1388 - t1787) * t1574 + (t1337 * t1389 - t1785) * t1577 + ((t1213 * t1777 - t1237 * t1330) * t1489 + t1697 * t1380) * t1403 + ((t1211 * t1779 - t1235 * t1328) * t1486 + t1698 * t1378) * t1400 + ((t1212 * t1781 - t1236 * t1326) * t1483 + t1699 * t1376) * t1397, t1376 * t1591 + t1378 * t1590 + t1380 * t1589 + (-t1380 * t1207 + (t1339 * t1744 + t1783) * t1795) * t1403 + (-t1378 * t1205 + (t1337 * t1745 + t1785) * t1799) * t1400 + (-t1376 * t1206 + (t1335 * t1746 + t1787) * t1803) * t1397, (t1700 * t1380 + (t1306 * t1330 + t1694 * t1777) * t1489) * t1403 + (t1701 * t1378 + (t1305 * t1328 + t1696 * t1779) * t1486) * t1400 + (t1702 * t1376 + (t1304 * t1326 + t1695 * t1781) * t1483) * t1397, (t1709 * t1380 + (t1330 * t1703 + t1706 * t1777) * t1489) * t1403 + (t1710 * t1378 + (t1328 * t1704 + t1707 * t1779) * t1486) * t1400 + (t1711 * t1376 + (t1326 * t1705 + t1708 * t1781) * t1483) * t1397, 0; 0, 0, 0, -t1454 * t1576 - t1455 * t1579 - t1456 * t1573, t1331 * t1727 + t1332 * t1725 + t1333 * t1723, t1520, t1510, t1235 * t1725 + t1236 * t1727 + t1237 * t1723, -t1520 * pkin(5) + (t1346 * t1733 + t1347 * t1734 + t1348 * t1732) * pkin(1), -t1510 * pkin(5) + (t1468 * t1576 + t1470 * t1579 + t1472 * t1573) * pkin(1), (-t1237 * t1411 + (t1213 + (-t1411 * t1476 + t1390) * t1348) * t1472) * t1489 + (-t1235 * t1410 + (t1211 + (-t1410 * t1475 + t1389) * t1346) * t1470) * t1486 + (-t1236 * t1409 + (t1212 + (-t1409 * t1474 + t1388) * t1347) * t1468) * t1483, (t1411 + t1416) * t1249 * t1723 + (t1410 + t1415) * t1248 * t1725 + (t1409 + t1414) * t1247 * t1727, (t1306 * t1411 + t1472 * t1694) * t1489 + (t1305 * t1410 + t1470 * t1696) * t1486 + (t1304 * t1409 + t1468 * t1695) * t1483, (t1411 * t1703 + t1472 * t1706) * t1489 + (t1410 * t1704 + t1470 * t1707) * t1486 + (t1409 * t1705 + t1468 * t1708) * t1483, 0;];
tau_reg  = t1;
