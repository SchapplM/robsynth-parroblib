% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRRRR1V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x18]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:34:42
% EndTime: 2020-08-07 03:35:00
% DurationCPUTime: 20.02s
% Computational Cost: add. (105495->712), mult. (210897->1440), div. (16155->16), fcn. (163608->36), ass. (0->608)
t1494 = qJ(2,3) + qJ(3,3);
t1440 = cos(qJ(1,3) + t1494);
t1441 = cos(qJ(1,3) - t1494);
t1422 = t1441 + t1440;
t1495 = qJ(2,2) + qJ(3,2);
t1442 = cos(qJ(1,2) + t1495);
t1443 = cos(qJ(1,2) - t1495);
t1423 = t1443 + t1442;
t1496 = qJ(2,1) + qJ(3,1);
t1444 = cos(qJ(1,1) + t1496);
t1445 = cos(qJ(1,1) - t1496);
t1424 = t1445 + t1444;
t1930 = pkin(2) ^ 2;
t1509 = cos(qJ(3,3));
t1485 = t1509 ^ 2;
t1512 = cos(qJ(3,2));
t1488 = t1512 ^ 2;
t1515 = cos(qJ(3,1));
t1491 = t1515 ^ 2;
t1902 = pkin(3) * t1509;
t1452 = pkin(2) + t1902;
t1497 = legFrame(3,2);
t1464 = sin(t1497);
t1467 = cos(t1497);
t1500 = sin(qJ(3,3));
t1501 = sin(qJ(2,3));
t1502 = sin(qJ(1,3));
t1510 = cos(qJ(2,3));
t1519 = xDP(2);
t1520 = xDP(1);
t1898 = pkin(3) * t1520;
t1764 = t1500 * t1898;
t1899 = pkin(3) * t1519;
t1765 = t1500 * t1899;
t1511 = cos(qJ(1,3));
t1518 = xDP(3);
t1795 = t1511 * t1518;
t1814 = t1502 * t1520;
t1815 = t1502 * t1519;
t1347 = ((-t1452 * t1814 + t1765) * t1467 + (t1452 * t1815 + t1764) * t1464 - t1452 * t1795) * t1510 + t1501 * ((t1452 * t1519 + t1502 * t1764) * t1467 + (t1452 * t1520 - t1502 * t1765) * t1464 + t1500 * pkin(3) * t1795);
t1777 = 0.2e1 * t1467;
t1796 = t1510 * t1520;
t1797 = t1510 * t1519;
t1380 = (-t1501 * t1814 - t1797) * t1777 + 0.2e1 * t1464 * (t1501 * t1815 - t1796);
t1381 = (-t1501 * t1519 + t1502 * t1796) * t1467 - t1464 * (t1501 * t1520 + t1502 * t1797);
t1356 = t1380 * t1500 + 0.2e1 * t1381 * t1509 + t1422 * t1518;
t1461 = sin(t1494);
t1476 = 0.1e1 / t1500;
t1525 = 0.1e1 / pkin(2);
t1395 = -t1502 * t1518 + (-t1464 * t1519 + t1467 * t1520) * t1511;
t1471 = t1510 * pkin(2);
t1455 = pkin(1) + t1471;
t1425 = pkin(3) * cos(t1494) + t1455;
t1417 = 0.1e1 / t1425 ^ 2;
t1862 = t1395 * t1417;
t1880 = t1356 * t1476;
t1671 = t1880 / 0.2e1;
t1350 = t1525 * t1671;
t1522 = 0.1e1 / pkin(3);
t1788 = t1522 * t1525;
t1688 = t1476 * t1788;
t1652 = t1347 * t1688;
t1344 = t1350 + t1652;
t1897 = t1344 * pkin(3);
t1907 = t1525 / 0.2e1;
t1302 = (t1461 * t1897 + (t1461 * t1347 * t1525 + (t1501 / 0.2e1 + (pkin(2) * t1501 + pkin(3) * t1461) * t1907) * t1356) * t1476) * t1862;
t1929 = -0.2e1 * t1302;
t1928 = 0.2e1 * t1302;
t1901 = pkin(3) * t1512;
t1453 = pkin(2) + t1901;
t1498 = legFrame(2,2);
t1465 = sin(t1498);
t1468 = cos(t1498);
t1503 = sin(qJ(3,2));
t1504 = sin(qJ(2,2));
t1505 = sin(qJ(1,2));
t1513 = cos(qJ(2,2));
t1760 = t1503 * t1898;
t1761 = t1503 * t1899;
t1514 = cos(qJ(1,2));
t1792 = t1514 * t1518;
t1806 = t1505 * t1520;
t1807 = t1505 * t1519;
t1348 = ((-t1453 * t1806 + t1761) * t1468 + (t1453 * t1807 + t1760) * t1465 - t1453 * t1792) * t1513 + t1504 * ((t1453 * t1519 + t1505 * t1760) * t1468 + (t1453 * t1520 - t1505 * t1761) * t1465 + t1503 * pkin(3) * t1792);
t1776 = 0.2e1 * t1468;
t1793 = t1513 * t1520;
t1794 = t1513 * t1519;
t1382 = (-t1504 * t1806 - t1794) * t1776 + 0.2e1 * t1465 * (t1504 * t1807 - t1793);
t1383 = (-t1504 * t1519 + t1505 * t1793) * t1468 - t1465 * (t1504 * t1520 + t1505 * t1794);
t1357 = t1382 * t1503 + 0.2e1 * t1383 * t1512 + t1423 * t1518;
t1462 = sin(t1495);
t1479 = 0.1e1 / t1503;
t1396 = -t1505 * t1518 + (-t1465 * t1519 + t1468 * t1520) * t1514;
t1473 = t1513 * pkin(2);
t1456 = pkin(1) + t1473;
t1426 = pkin(3) * cos(t1495) + t1456;
t1419 = 0.1e1 / t1426 ^ 2;
t1860 = t1396 * t1419;
t1878 = t1357 * t1479;
t1669 = t1878 / 0.2e1;
t1351 = t1525 * t1669;
t1686 = t1479 * t1788;
t1651 = t1348 * t1686;
t1345 = t1351 + t1651;
t1896 = t1345 * pkin(3);
t1303 = (t1462 * t1896 + (t1462 * t1348 * t1525 + (t1504 / 0.2e1 + (pkin(2) * t1504 + pkin(3) * t1462) * t1907) * t1357) * t1479) * t1860;
t1927 = -0.2e1 * t1303;
t1926 = 0.2e1 * t1303;
t1900 = pkin(3) * t1515;
t1454 = pkin(2) + t1900;
t1499 = legFrame(1,2);
t1466 = sin(t1499);
t1469 = cos(t1499);
t1506 = sin(qJ(3,1));
t1507 = sin(qJ(2,1));
t1508 = sin(qJ(1,1));
t1516 = cos(qJ(2,1));
t1756 = t1506 * t1898;
t1757 = t1506 * t1899;
t1517 = cos(qJ(1,1));
t1789 = t1517 * t1518;
t1798 = t1508 * t1520;
t1799 = t1508 * t1519;
t1349 = ((-t1454 * t1798 + t1757) * t1469 + (t1454 * t1799 + t1756) * t1466 - t1454 * t1789) * t1516 + t1507 * ((t1454 * t1519 + t1508 * t1756) * t1469 + (t1454 * t1520 - t1508 * t1757) * t1466 + t1506 * pkin(3) * t1789);
t1775 = 0.2e1 * t1469;
t1790 = t1516 * t1520;
t1791 = t1516 * t1519;
t1384 = (-t1507 * t1798 - t1791) * t1775 + 0.2e1 * t1466 * (t1507 * t1799 - t1790);
t1385 = (-t1507 * t1519 + t1508 * t1790) * t1469 - t1466 * (t1507 * t1520 + t1508 * t1791);
t1358 = t1384 * t1506 + 0.2e1 * t1385 * t1515 + t1424 * t1518;
t1463 = sin(t1496);
t1482 = 0.1e1 / t1506;
t1397 = -t1508 * t1518 + (-t1466 * t1519 + t1469 * t1520) * t1517;
t1475 = t1516 * pkin(2);
t1457 = pkin(1) + t1475;
t1427 = pkin(3) * cos(t1496) + t1457;
t1421 = 0.1e1 / t1427 ^ 2;
t1858 = t1397 * t1421;
t1876 = t1358 * t1482;
t1667 = t1876 / 0.2e1;
t1352 = t1525 * t1667;
t1684 = t1482 * t1788;
t1650 = t1349 * t1684;
t1346 = t1352 + t1650;
t1895 = t1346 * pkin(3);
t1304 = (t1463 * t1895 + (t1463 * t1349 * t1525 + (t1507 / 0.2e1 + (pkin(2) * t1507 + pkin(3) * t1463) * t1907) * t1358) * t1482) * t1858;
t1925 = -0.2e1 * t1304;
t1924 = 0.2e1 * t1304;
t1470 = t1509 * pkin(2);
t1923 = (t1470 + pkin(3)) * t1522;
t1472 = t1512 * pkin(2);
t1922 = (t1472 + pkin(3)) * t1522;
t1474 = t1515 * pkin(2);
t1921 = (t1474 + pkin(3)) * t1522;
t1486 = t1510 ^ 2;
t1784 = t1485 - 0.1e1 / 0.2e1;
t1817 = t1501 * t1510;
t1392 = t1784 * t1817 + (t1486 - 0.1e1 / 0.2e1) * t1509 * t1500;
t1920 = 0.4e1 * t1392;
t1489 = t1513 ^ 2;
t1783 = t1488 - 0.1e1 / 0.2e1;
t1809 = t1504 * t1513;
t1393 = t1783 * t1809 + (t1489 - 0.1e1 / 0.2e1) * t1512 * t1503;
t1919 = 0.4e1 * t1393;
t1492 = t1516 ^ 2;
t1782 = t1491 - 0.1e1 / 0.2e1;
t1801 = t1507 * t1516;
t1394 = t1782 * t1801 + (t1492 - 0.1e1 / 0.2e1) * t1515 * t1506;
t1918 = 0.4e1 * t1394;
t1917 = 0.2e1 * t1471;
t1916 = 0.2e1 * t1473;
t1915 = 0.2e1 * t1475;
t1914 = -0.2e1 * t1486;
t1913 = -0.2e1 * t1489;
t1912 = -0.2e1 * t1492;
t1521 = pkin(3) ^ 2;
t1911 = -0.2e1 * t1521;
t1910 = t1422 / 0.2e1;
t1909 = t1423 / 0.2e1;
t1908 = t1424 / 0.2e1;
t1906 = t1525 / 0.4e1;
t1458 = 0.2e1 * t1486 - 0.1e1;
t1459 = 0.2e1 * t1489 - 0.1e1;
t1460 = 0.2e1 * t1492 - 0.1e1;
t1905 = pkin(1) * t1510;
t1904 = pkin(1) * t1513;
t1903 = pkin(1) * t1516;
t1894 = t1485 * pkin(3);
t1893 = t1488 * pkin(3);
t1892 = t1491 * pkin(3);
t1891 = -0.2e1 * pkin(2) * pkin(3);
t1770 = t1509 * t1897;
t1326 = t1671 + t1770;
t1335 = -0.2e1 * t1897;
t1820 = t1500 * t1510;
t1766 = pkin(3) * t1820;
t1401 = t1452 * t1501 + t1766;
t1428 = -pkin(3) + t1470 + 0.2e1 * t1894;
t1487 = t1511 ^ 2;
t1830 = t1476 * t1509;
t1680 = t1356 * t1830;
t1683 = t1500 * t1817;
t1416 = 0.1e1 / t1425;
t1853 = t1416 * t1511;
t1712 = t1395 * t1853;
t1821 = t1500 * t1501;
t1767 = pkin(3) * t1821;
t1780 = 0.2e1 * t1902;
t1854 = t1416 * t1502;
t1890 = (((0.4e1 * t1344 * t1894 + t1335 + t1680) * t1486 - 0.2e1 * (t1671 + 0.2e1 * t1770) * t1683 - t1335 + t1335 * t1485) * t1487 - 0.2e1 * t1502 * (t1428 * t1817 + ((pkin(2) + t1780) * t1486 - t1452) * t1500) * t1712 - t1680 + t1335 + t1422 * ((-t1326 * t1510 + t1344 * t1767) * t1511 + t1401 * t1395 * t1854)) * t1356;
t1769 = t1512 * t1896;
t1327 = t1669 + t1769;
t1336 = -0.2e1 * t1896;
t1812 = t1503 * t1513;
t1762 = pkin(3) * t1812;
t1402 = t1453 * t1504 + t1762;
t1429 = -pkin(3) + t1472 + 0.2e1 * t1893;
t1490 = t1514 ^ 2;
t1828 = t1479 * t1512;
t1678 = t1357 * t1828;
t1682 = t1503 * t1809;
t1418 = 0.1e1 / t1426;
t1850 = t1418 * t1514;
t1710 = t1396 * t1850;
t1813 = t1503 * t1504;
t1763 = pkin(3) * t1813;
t1779 = 0.2e1 * t1901;
t1851 = t1418 * t1505;
t1889 = (((0.4e1 * t1345 * t1893 + t1336 + t1678) * t1489 - 0.2e1 * (t1669 + 0.2e1 * t1769) * t1682 - t1336 + t1336 * t1488) * t1490 - 0.2e1 * t1505 * (t1429 * t1809 + ((pkin(2) + t1779) * t1489 - t1453) * t1503) * t1710 - t1678 + t1336 + t1423 * ((-t1327 * t1513 + t1345 * t1763) * t1514 + t1402 * t1396 * t1851)) * t1357;
t1768 = t1515 * t1895;
t1328 = t1667 + t1768;
t1337 = -0.2e1 * t1895;
t1804 = t1506 * t1516;
t1758 = pkin(3) * t1804;
t1403 = t1454 * t1507 + t1758;
t1430 = -pkin(3) + t1474 + 0.2e1 * t1892;
t1493 = t1517 ^ 2;
t1826 = t1482 * t1515;
t1676 = t1358 * t1826;
t1681 = t1506 * t1801;
t1420 = 0.1e1 / t1427;
t1847 = t1420 * t1517;
t1708 = t1397 * t1847;
t1805 = t1506 * t1507;
t1759 = pkin(3) * t1805;
t1778 = 0.2e1 * t1900;
t1848 = t1420 * t1508;
t1888 = (((0.4e1 * t1346 * t1892 + t1337 + t1676) * t1492 - 0.2e1 * (t1667 + 0.2e1 * t1768) * t1681 - t1337 + t1337 * t1491) * t1493 - 0.2e1 * t1508 * (t1430 * t1801 + ((pkin(2) + t1778) * t1492 - t1454) * t1506) * t1708 - t1676 + t1337 + t1424 * ((-t1328 * t1516 + t1346 * t1759) * t1517 + t1403 * t1397 * t1848)) * t1358;
t1887 = t1302 * t1416;
t1886 = t1302 * t1476;
t1885 = t1303 * t1418;
t1884 = t1303 * t1479;
t1883 = t1304 * t1420;
t1882 = t1304 * t1482;
t1881 = t1356 * t1395;
t1879 = t1357 * t1396;
t1877 = t1358 * t1397;
t1410 = -t1509 * t1510 + t1821;
t1818 = t1501 * t1509;
t1413 = t1818 + t1820;
t1842 = t1464 * t1502;
t1368 = t1410 * t1842 - t1413 * t1467;
t1875 = t1368 * t1476;
t1836 = t1467 * t1502;
t1369 = -t1410 * t1836 - t1413 * t1464;
t1874 = t1369 * t1476;
t1411 = -t1512 * t1513 + t1813;
t1810 = t1504 * t1512;
t1414 = t1810 + t1812;
t1840 = t1465 * t1505;
t1370 = t1411 * t1840 - t1414 * t1468;
t1873 = t1370 * t1479;
t1834 = t1468 * t1505;
t1371 = -t1411 * t1834 - t1414 * t1465;
t1872 = t1371 * t1479;
t1412 = -t1515 * t1516 + t1805;
t1802 = t1507 * t1515;
t1415 = t1802 + t1804;
t1838 = t1466 * t1508;
t1372 = t1412 * t1838 - t1415 * t1469;
t1871 = t1372 * t1482;
t1832 = t1469 * t1508;
t1373 = -t1412 * t1832 - t1415 * t1466;
t1870 = t1373 * t1482;
t1431 = t1452 * t1510;
t1404 = t1431 - t1767;
t1374 = t1401 * t1464 - t1404 * t1836;
t1869 = t1374 * t1476;
t1432 = t1453 * t1513;
t1405 = t1432 - t1763;
t1375 = t1402 * t1465 - t1405 * t1834;
t1868 = t1375 * t1479;
t1433 = t1454 * t1516;
t1406 = t1433 - t1759;
t1376 = t1403 * t1466 - t1406 * t1832;
t1867 = t1376 * t1482;
t1377 = t1401 * t1467 + t1404 * t1842;
t1866 = t1377 * t1476;
t1378 = t1402 * t1468 + t1405 * t1840;
t1865 = t1378 * t1479;
t1379 = t1403 * t1469 + t1406 * t1838;
t1864 = t1379 * t1482;
t1389 = t1395 ^ 2;
t1386 = t1389 * t1417;
t1390 = t1396 ^ 2;
t1387 = t1390 * t1419;
t1391 = t1397 ^ 2;
t1388 = t1391 * t1421;
t1863 = t1395 * t1416;
t1861 = t1396 * t1418;
t1859 = t1397 * t1420;
t1857 = t1404 * t1511;
t1856 = t1405 * t1514;
t1855 = t1406 * t1517;
t1852 = t1417 * t1476;
t1849 = t1419 * t1479;
t1846 = t1421 * t1482;
t1845 = t1422 * t1476;
t1844 = t1423 * t1479;
t1843 = t1424 * t1482;
t1841 = t1464 * t1511;
t1839 = t1465 * t1514;
t1837 = t1466 * t1517;
t1835 = t1467 * t1511;
t1833 = t1468 * t1514;
t1831 = t1469 * t1517;
t1477 = 0.1e1 / t1500 ^ 2;
t1526 = 0.1e1 / t1930;
t1829 = t1477 * t1526;
t1480 = 0.1e1 / t1503 ^ 2;
t1827 = t1480 * t1526;
t1483 = 0.1e1 / t1506 ^ 2;
t1825 = t1483 * t1526;
t1824 = t1486 * t1500;
t1823 = t1489 * t1503;
t1822 = t1492 * t1506;
t1819 = t1501 * t1502;
t1816 = t1502 * t1510;
t1811 = t1504 * t1505;
t1808 = t1505 * t1513;
t1803 = t1507 * t1508;
t1800 = t1508 * t1516;
t1338 = t1350 + t1652 / 0.2e1;
t1434 = pkin(1) - 0.2e1 * t1767;
t1755 = pkin(1) * t1821;
t1577 = -pkin(3) + t1755 + t1894;
t1662 = -t1522 * t1526 / 0.2e1;
t1781 = t1521 - t1930;
t1787 = (t1344 * t1521 + (t1338 * t1780 + t1671) * pkin(2)) * t1477 * t1356 * t1662 + ((t1485 * t1911 + t1509 * t1891 + t1781) * t1486 - t1434 * t1431 + pkin(3) * t1577) * t1688 * t1386;
t1339 = t1351 + t1651 / 0.2e1;
t1435 = pkin(1) - 0.2e1 * t1763;
t1754 = pkin(1) * t1813;
t1576 = -pkin(3) + t1754 + t1893;
t1786 = (t1345 * t1521 + (t1339 * t1779 + t1669) * pkin(2)) * t1480 * t1357 * t1662 + ((t1488 * t1911 + t1512 * t1891 + t1781) * t1489 - t1435 * t1432 + pkin(3) * t1576) * t1686 * t1387;
t1340 = t1352 + t1650 / 0.2e1;
t1436 = pkin(1) - 0.2e1 * t1759;
t1753 = pkin(1) * t1805;
t1575 = -pkin(3) + t1753 + t1892;
t1785 = (t1346 * t1521 + (t1340 * t1778 + t1667) * pkin(2)) * t1483 * t1358 * t1662 + ((t1491 * t1911 + t1515 * t1891 + t1781) * t1492 - t1436 * t1433 + pkin(3) * t1575) * t1684 * t1388;
t1774 = 0.2e1 * t1522;
t1773 = pkin(2) * t1821;
t1772 = pkin(2) * t1813;
t1771 = pkin(2) * t1805;
t1725 = t1347 * t1829;
t1323 = t1344 * t1725;
t1670 = t1880 / 0.4e1;
t1604 = (-0.4e1 * t1502 * ((t1670 + t1770) * t1824 + (t1509 * t1670 + t1784 * t1897) * t1817 - t1500 * t1326 / 0.2e1) * t1511 + (0.2e1 * t1487 - 0.2e1) * (t1428 * t1486 + (t1434 * t1509 - t1773) * t1510 - t1577) * t1863 + t1422 * ((t1326 * t1501 + t1344 * t1766) * t1502 - (pkin(1) + t1404) * t1712)) * t1476 * t1525 * t1863;
t1281 = -t1604 / 0.2e1 - t1829 * t1890 / 0.4e1 + t1323;
t1752 = t1281 * t1830;
t1724 = t1348 * t1827;
t1324 = t1345 * t1724;
t1668 = t1878 / 0.4e1;
t1603 = (-0.4e1 * t1505 * ((t1668 + t1769) * t1823 + (t1512 * t1668 + t1783 * t1896) * t1809 - t1503 * t1327 / 0.2e1) * t1514 + (0.2e1 * t1490 - 0.2e1) * (t1429 * t1489 + (t1435 * t1512 - t1772) * t1513 - t1576) * t1861 + t1423 * ((t1327 * t1504 + t1345 * t1762) * t1505 - (pkin(1) + t1405) * t1710)) * t1479 * t1525 * t1861;
t1282 = -t1603 / 0.2e1 - t1827 * t1889 / 0.4e1 + t1324;
t1751 = t1282 * t1828;
t1723 = t1349 * t1825;
t1325 = t1346 * t1723;
t1666 = t1876 / 0.4e1;
t1602 = (-0.4e1 * t1508 * ((t1666 + t1768) * t1822 + (t1515 * t1666 + t1782 * t1895) * t1801 - t1506 * t1328 / 0.2e1) * t1517 + (0.2e1 * t1493 - 0.2e1) * (t1430 * t1492 + (t1436 * t1515 - t1771) * t1516 - t1575) * t1859 + t1424 * ((t1328 * t1507 + t1346 * t1758) * t1508 - (pkin(1) + t1406) * t1708)) * t1482 * t1525 * t1859;
t1283 = -t1602 / 0.2e1 - t1825 * t1888 / 0.4e1 + t1325;
t1750 = t1283 * t1826;
t1749 = t1410 * t1886;
t1748 = t1413 * t1886;
t1747 = t1455 * t1887;
t1746 = t1302 * t1854;
t1745 = t1302 * t1853;
t1744 = t1501 * t1886;
t1743 = t1510 * t1886;
t1742 = t1411 * t1884;
t1741 = t1414 * t1884;
t1740 = t1456 * t1885;
t1739 = t1303 * t1851;
t1738 = t1303 * t1850;
t1737 = t1504 * t1884;
t1736 = t1513 * t1884;
t1735 = t1412 * t1882;
t1734 = t1415 * t1882;
t1733 = t1457 * t1883;
t1732 = t1304 * t1848;
t1731 = t1304 * t1847;
t1730 = t1507 * t1882;
t1729 = t1516 * t1882;
t1728 = t1344 * t1862;
t1727 = t1345 * t1860;
t1726 = t1346 * t1858;
t1353 = t1356 ^ 2;
t1722 = t1353 * t1416 * t1477;
t1354 = t1357 ^ 2;
t1721 = t1354 * t1418 * t1480;
t1355 = t1358 ^ 2;
t1720 = t1355 * t1420 * t1483;
t1719 = t1511 * t1881;
t1718 = t1514 * t1879;
t1717 = t1517 * t1877;
t1716 = t1389 * t1852;
t1715 = t1390 * t1849;
t1714 = t1391 * t1846;
t1707 = t1476 * t1857;
t1706 = t1479 * t1856;
t1705 = t1482 * t1855;
t1704 = t1416 * t1835;
t1703 = t1416 * t1819;
t1702 = t1501 * t1853;
t1701 = t1416 * t1816;
t1700 = t1510 * t1853;
t1699 = t1418 * t1833;
t1698 = t1418 * t1811;
t1697 = t1504 * t1850;
t1696 = t1418 * t1808;
t1695 = t1513 * t1850;
t1694 = t1420 * t1831;
t1693 = t1420 * t1803;
t1692 = t1507 * t1847;
t1691 = t1420 * t1800;
t1690 = t1516 * t1847;
t1278 = -t1604 + 0.2e1 * t1323 + (-t1890 / 0.2e1 - t1344 * t1347 * t1923) * t1829 + t1787;
t1689 = t1278 * t1830;
t1279 = -t1603 + 0.2e1 * t1324 + (-t1889 / 0.2e1 - t1345 * t1348 * t1922) * t1827 + t1786;
t1687 = t1279 * t1828;
t1280 = -t1602 + 0.2e1 * t1325 + (-t1888 / 0.2e1 - t1346 * t1349 * t1921) * t1825 + t1785;
t1685 = t1280 * t1826;
t1679 = t1281 * t1857;
t1677 = t1282 * t1856;
t1675 = t1283 * t1855;
t1674 = t1353 * t1906;
t1673 = t1354 * t1906;
t1672 = t1355 * t1906;
t1665 = t1845 / 0.2e1;
t1664 = t1844 / 0.2e1;
t1663 = t1843 / 0.2e1;
t1661 = t1511 * t1777;
t1660 = t1514 * t1776;
t1659 = t1517 * t1775;
t1658 = t1302 * t1707;
t1478 = t1501 ^ 2;
t1657 = t1478 * t1745;
t1656 = t1303 * t1706;
t1481 = t1504 ^ 2;
t1655 = t1481 * t1738;
t1654 = t1304 * t1705;
t1484 = t1507 ^ 2;
t1653 = t1484 * t1731;
t1649 = t1356 * t1458 * t1862;
t1648 = t1357 * t1459 * t1860;
t1647 = t1358 * t1460 * t1858;
t1646 = t1392 * t1716;
t1398 = -0.4e1 * t1509 * t1683 + (0.4e1 * t1486 - 0.2e1) * t1485 - t1458;
t1645 = t1398 * t1716;
t1644 = t1422 * t1716;
t1643 = (t1917 + pkin(1)) * t1501 * t1386;
t1642 = t1393 * t1715;
t1399 = -0.4e1 * t1512 * t1682 + (0.4e1 * t1489 - 0.2e1) * t1488 - t1459;
t1641 = t1399 * t1715;
t1640 = t1423 * t1715;
t1639 = (t1916 + pkin(1)) * t1504 * t1387;
t1638 = t1394 * t1714;
t1400 = -0.4e1 * t1515 * t1681 + (0.4e1 * t1492 - 0.2e1) * t1491 - t1460;
t1637 = t1400 * t1714;
t1636 = (t1915 + pkin(1)) * t1507 * t1388;
t1635 = t1464 * t1702;
t1634 = t1464 * t1700;
t1633 = t1467 * t1702;
t1632 = t1467 * t1700;
t1631 = t1817 * t1852;
t1630 = t1465 * t1697;
t1629 = t1465 * t1695;
t1628 = t1468 * t1697;
t1627 = t1468 * t1695;
t1626 = t1809 * t1849;
t1625 = t1466 * t1692;
t1624 = t1466 * t1690;
t1623 = t1469 * t1692;
t1622 = t1469 * t1690;
t1621 = t1801 * t1846;
t1620 = t1424 * t1714;
t1619 = -t1302 * t1845 / 0.2e1;
t1618 = -t1303 * t1844 / 0.2e1;
t1617 = -t1304 * t1843 / 0.2e1;
t1616 = -t1722 / 0.4e1;
t1615 = t1722 / 0.4e1;
t1614 = -t1721 / 0.4e1;
t1613 = t1721 / 0.4e1;
t1612 = -t1720 / 0.4e1;
t1611 = t1720 / 0.4e1;
t1610 = -0.2e1 * t1646;
t1609 = -0.2e1 * t1642;
t1608 = -0.2e1 * t1638;
t1607 = t1701 * t1929;
t1606 = t1696 * t1927;
t1605 = t1691 * t1925;
t1362 = (t1914 + 0.1e1) * t1386;
t1363 = (t1913 + 0.1e1) * t1387;
t1364 = (t1912 + 0.1e1) * t1388;
t1601 = t1511 * t1649;
t1600 = t1514 * t1648;
t1599 = t1517 * t1647;
t1598 = t1707 * t1386;
t1597 = t1706 * t1387;
t1596 = t1705 * t1388;
t1595 = t1511 * t1616;
t1594 = t1511 * t1615;
t1593 = t1514 * t1614;
t1592 = t1514 * t1613;
t1591 = t1517 * t1612;
t1590 = t1517 * t1611;
t1275 = -t1323 * t1923 + t1281 + t1787;
t1341 = t1344 ^ 2;
t1589 = t1275 * t1413 - t1341 * t1410;
t1588 = t1275 * t1410 + t1341 * t1413;
t1276 = -t1324 * t1922 + t1282 + t1786;
t1342 = t1345 ^ 2;
t1587 = t1276 * t1414 - t1342 * t1411;
t1586 = t1276 * t1411 + t1342 * t1414;
t1277 = -t1325 * t1921 + t1283 + t1785;
t1343 = t1346 ^ 2;
t1585 = t1277 * t1415 - t1343 * t1412;
t1584 = t1277 * t1412 + t1343 * t1415;
t1583 = t1634 * t1929;
t1582 = t1629 * t1927;
t1581 = t1624 * t1925;
t1580 = t1632 * t1928;
t1579 = t1627 * t1926;
t1578 = t1622 * t1924;
t1574 = pkin(2) * t1338 * t1914 - t1344 * t1905;
t1573 = pkin(2) * t1339 * t1913 - t1345 * t1904;
t1572 = pkin(2) * t1340 * t1912 - t1346 * t1903;
t1571 = t1416 * t1589;
t1570 = t1416 * t1588;
t1569 = t1418 * t1587;
t1568 = t1418 * t1586;
t1567 = t1420 * t1585;
t1566 = t1420 * t1584;
t1565 = t1389 * t1910 + t1502 * t1881;
t1564 = t1390 * t1909 + t1505 * t1879;
t1563 = t1391 * t1908 + t1508 * t1877;
t1562 = (pkin(1) * t1344 + t1338 * t1917) * t1501;
t1561 = (pkin(1) * t1345 + t1339 * t1916) * t1504;
t1560 = (pkin(1) * t1346 + t1340 * t1915) * t1507;
t1559 = t1368 * t1389 + t1464 * t1719;
t1558 = -t1369 * t1389 + t1467 * t1719;
t1557 = t1370 * t1390 + t1465 * t1718;
t1556 = -t1371 * t1390 + t1468 * t1718;
t1555 = t1372 * t1391 + t1466 * t1717;
t1554 = -t1373 * t1391 + t1469 * t1717;
t1553 = -(t1500 * t1562 + (t1671 + t1574) * t1509) * t1862 + t1413 * t1747;
t1552 = -(t1503 * t1561 + (t1669 + t1573) * t1512) * t1860 + t1414 * t1740;
t1551 = -(t1506 * t1560 + (t1667 + t1572) * t1515) * t1858 + t1415 * t1733;
t1550 = ((-t1443 / 0.2e1 - t1442 / 0.2e1) * t1518 + (-t1383 + t1561) * t1512 + (-t1382 / 0.2e1 - t1573) * t1503) * t1860 + t1411 * t1740;
t1549 = ((-t1441 / 0.2e1 - t1440 / 0.2e1) * t1518 + (-t1381 + t1562) * t1509 + (-t1380 / 0.2e1 - t1574) * t1500) * t1862 + t1410 * t1747;
t1548 = ((-t1445 / 0.2e1 - t1444 / 0.2e1) * t1518 + (-t1385 + t1560) * t1515 + (-t1384 / 0.2e1 - t1572) * t1506) * t1858 + t1412 * t1733;
t1547 = 0.2e1 * t1398 * t1728 + t1887 * t1920;
t1546 = t1413 ^ 2 * t1887 + t1728 * t1920;
t1545 = 0.2e1 * t1399 * t1727 + t1885 * t1919;
t1544 = t1414 ^ 2 * t1885 + t1727 * t1919;
t1543 = 0.2e1 * t1400 * t1726 + t1883 * t1918;
t1542 = t1415 ^ 2 * t1883 + t1726 * t1918;
t1541 = t1565 * t1852;
t1540 = t1564 * t1849;
t1539 = t1563 * t1846;
t1538 = t1559 * t1852;
t1537 = t1558 * t1852;
t1536 = t1557 * t1849;
t1535 = t1556 * t1849;
t1534 = t1555 * t1846;
t1533 = t1554 * t1846;
t1532 = 0.2e1 * t1553;
t1531 = 0.2e1 * t1552;
t1530 = 0.2e1 * t1551;
t1529 = 0.2e1 * t1550;
t1528 = 0.2e1 * t1549;
t1527 = 0.2e1 * t1548;
t1523 = 0.1e1 / pkin(3) ^ 2;
t1334 = t1482 * t1672 + (pkin(2) * t1822 + (pkin(1) * t1506 + pkin(2) * t1802) * t1516 + pkin(1) * t1802) * t1388;
t1333 = t1479 * t1673 + (pkin(2) * t1823 + (pkin(1) * t1503 + pkin(2) * t1810) * t1513 + pkin(1) * t1810) * t1387;
t1332 = t1476 * t1674 + (pkin(2) * t1824 + (pkin(1) * t1500 + pkin(2) * t1818) * t1510 + pkin(1) * t1818) * t1386;
t1331 = t1515 * t1483 * t1672 + (t1492 * t1474 + (pkin(1) * t1515 - t1771) * t1516 - t1753) * t1388;
t1330 = t1512 * t1480 * t1673 + (t1489 * t1472 + (pkin(1) * t1512 - t1772) * t1513 - t1754) * t1387;
t1329 = t1509 * t1477 * t1674 + (t1486 * t1470 + (pkin(1) * t1509 - t1773) * t1510 - t1755) * t1386;
t1322 = t1388 * t1903 + ((-t1349 * t1523 - t1358 * t1522) * t1723 - t1364) * pkin(2);
t1321 = t1387 * t1904 + ((-t1348 * t1523 - t1357 * t1522) * t1724 - t1363) * pkin(2);
t1320 = t1386 * t1905 + ((-t1347 * t1523 - t1356 * t1522) * t1725 - t1362) * pkin(2);
t1319 = t1322 * t1515 - t1506 * t1636;
t1318 = t1322 * t1506 + t1515 * t1636;
t1317 = t1321 * t1512 - t1503 * t1639;
t1316 = t1321 * t1503 + t1512 * t1639;
t1315 = t1320 * t1509 - t1500 * t1643;
t1314 = t1320 * t1500 + t1509 * t1643;
t1 = [t1302 * t1704 + t1303 * t1699 + t1304 * t1694, 0, 0, t1467 * t1657 + t1468 * t1655 + t1469 * t1653 + (t1554 * t1621 + t1556 * t1626 + t1558 * t1631) * t1525, t1501 * t1580 + t1504 * t1579 + t1507 * t1578 + ((t1364 * t1373 + t1469 * t1599) * t1482 + (t1363 * t1371 + t1468 * t1600) * t1479 + (t1362 * t1369 + t1467 * t1601) * t1476) * t1525, -t1281 * t1633 - t1282 * t1628 - t1283 * t1623 + (t1467 * t1510 * t1595 + t1468 * t1513 * t1593 + t1469 * t1516 * t1591) * t1526 + (-t1369 * t1744 - t1371 * t1737 - t1373 * t1730) * t1525, -t1281 * t1632 - t1282 * t1627 - t1283 * t1622 + (t1467 * t1501 * t1594 + t1468 * t1504 * t1592 + t1469 * t1507 * t1590) * t1526 + (-t1369 * t1743 - t1371 * t1736 - t1373 * t1729) * t1525, (t1281 * t1874 + t1282 * t1872 + t1283 * t1870) * t1525, (t1580 + t1579 + t1578 + (-t1501 * t1537 - t1504 * t1535 - t1507 * t1533) * t1525) * pkin(1), (t1633 * t1929 + t1628 * t1927 + t1623 * t1925 + (-t1510 * t1537 - t1513 * t1535 - t1516 * t1533) * t1525) * pkin(1), t1542 * t1831 + t1544 * t1833 + t1546 * t1835 + (t1369 * t1610 + t1371 * t1609 + t1373 * t1608 + (-t1374 * t1646 - t1375 * t1642 - t1376 * t1638) * t1774) * t1525, t1543 * t1831 + t1545 * t1833 + t1547 * t1835 + (-t1369 * t1645 - t1371 * t1641 - t1373 * t1637 + (-t1374 * t1645 - t1375 * t1641 - t1376 * t1637) * t1522) * t1525, -t1585 * t1694 - t1587 * t1699 - t1589 * t1704 + (-t1369 * t1748 - t1371 * t1741 - t1373 * t1734 + (-t1374 * t1748 - t1375 * t1741 - t1376 * t1734) * t1522) * t1525, t1584 * t1694 + t1586 * t1699 + t1588 * t1704 + (t1369 * t1749 + t1371 * t1742 + t1373 * t1735 + (t1374 * t1749 + t1375 * t1742 + t1376 * t1735) * t1522) * t1525, (t1275 * t1874 + t1276 * t1872 + t1277 * t1870 + (t1275 * t1869 + t1276 * t1868 + t1277 * t1867) * t1522) * t1525, t1369 * t1689 + t1371 * t1687 + t1373 * t1685 - t1548 * t1659 - t1550 * t1660 - t1549 * t1661 + (t1314 * t1874 + t1316 * t1872 + t1318 * t1870) * t1525 + (t1374 * t1752 + t1375 * t1751 + t1376 * t1750 + (t1332 * t1869 + t1333 * t1868 + t1334 * t1867) * t1525) * t1522, -t1369 * t1278 - t1371 * t1279 - t1373 * t1280 - t1551 * t1659 - t1552 * t1660 - t1553 * t1661 + (t1315 * t1874 + t1317 * t1872 + t1319 * t1870) * t1525 + (-t1374 * t1281 - t1375 * t1282 - t1376 * t1283 + (t1329 * t1869 + t1330 * t1868 + t1331 * t1867) * t1525) * t1522, 0; -t1464 * t1745 - t1465 * t1738 - t1466 * t1731, 0, 0, -t1464 * t1657 - t1465 * t1655 - t1466 * t1653 + (-t1555 * t1621 - t1557 * t1626 - t1559 * t1631) * t1525, t1501 * t1583 + t1504 * t1582 + t1507 * t1581 + ((t1364 * t1372 - t1466 * t1599) * t1482 + (t1363 * t1370 - t1465 * t1600) * t1479 + (t1362 * t1368 - t1464 * t1601) * t1476) * t1525, t1281 * t1635 + t1282 * t1630 + t1283 * t1625 + (t1464 * t1510 * t1594 + t1465 * t1513 * t1592 + t1466 * t1516 * t1590) * t1526 + (-t1368 * t1744 - t1370 * t1737 - t1372 * t1730) * t1525, t1281 * t1634 + t1282 * t1629 + t1283 * t1624 + (t1464 * t1501 * t1595 + t1465 * t1504 * t1593 + t1466 * t1507 * t1591) * t1526 + (-t1368 * t1743 - t1370 * t1736 - t1372 * t1729) * t1525, (t1281 * t1875 + t1282 * t1873 + t1283 * t1871) * t1525, (t1583 + t1582 + t1581 + (t1501 * t1538 + t1504 * t1536 + t1507 * t1534) * t1525) * pkin(1), (t1635 * t1928 + t1630 * t1926 + t1625 * t1924 + (t1510 * t1538 + t1513 * t1536 + t1516 * t1534) * t1525) * pkin(1), -t1542 * t1837 - t1544 * t1839 - t1546 * t1841 + (t1368 * t1610 + t1370 * t1609 + t1372 * t1608 + (-t1377 * t1646 - t1378 * t1642 - t1379 * t1638) * t1774) * t1525, -t1543 * t1837 - t1545 * t1839 - t1547 * t1841 + (-t1368 * t1645 - t1370 * t1641 - t1372 * t1637 + (-t1377 * t1645 - t1378 * t1641 - t1379 * t1637) * t1522) * t1525, t1567 * t1837 + t1569 * t1839 + t1571 * t1841 + (-t1368 * t1748 - t1370 * t1741 - t1372 * t1734 + (-t1377 * t1748 - t1378 * t1741 - t1379 * t1734) * t1522) * t1525, -t1566 * t1837 - t1568 * t1839 - t1570 * t1841 + (t1368 * t1749 + t1370 * t1742 + t1372 * t1735 + (t1377 * t1749 + t1378 * t1742 + t1379 * t1735) * t1522) * t1525, (t1275 * t1875 + t1276 * t1873 + t1277 * t1871 + (t1275 * t1866 + t1276 * t1865 + t1277 * t1864) * t1522) * t1525, t1368 * t1689 + t1370 * t1687 + t1372 * t1685 + t1527 * t1837 + t1529 * t1839 + t1528 * t1841 + (t1314 * t1875 + t1316 * t1873 + t1318 * t1871) * t1525 + (t1377 * t1752 + t1378 * t1751 + t1379 * t1750 + (t1332 * t1866 + t1333 * t1865 + t1334 * t1864) * t1525) * t1522, -t1368 * t1278 - t1370 * t1279 - t1372 * t1280 + t1530 * t1837 + t1531 * t1839 + t1532 * t1841 + (t1315 * t1875 + t1317 * t1873 + t1319 * t1871) * t1525 + (-t1377 * t1281 - t1378 * t1282 - t1379 * t1283 + (t1329 * t1866 + t1330 * t1865 + t1331 * t1864) * t1525) * t1522, 0; -t1732 - t1739 - t1746, 0, 0, -t1478 * t1746 - t1481 * t1739 - t1484 * t1732 + (-t1563 * t1621 - t1564 * t1626 - t1565 * t1631) * t1525, t1501 * t1607 + t1504 * t1606 + t1507 * t1605 + ((t1364 * t1908 - t1508 * t1647) * t1482 + (t1363 * t1909 - t1505 * t1648) * t1479 + (t1362 * t1910 - t1502 * t1649) * t1476) * t1525, t1281 * t1703 + t1282 * t1698 + t1283 * t1693 + (t1611 * t1800 + t1613 * t1808 + t1615 * t1816) * t1526 + (t1501 * t1619 + t1504 * t1618 + t1507 * t1617) * t1525, t1281 * t1701 + t1282 * t1696 + t1283 * t1691 + (t1612 * t1803 + t1614 * t1811 + t1616 * t1819) * t1526 + (t1510 * t1619 + t1513 * t1618 + t1516 * t1617) * t1525, (t1281 * t1665 + t1282 * t1664 + t1283 * t1663) * t1525, (t1607 + t1606 + t1605 + (t1501 * t1541 + t1504 * t1540 + t1507 * t1539) * t1525) * pkin(1), (t1703 * t1928 + t1698 * t1926 + t1693 * t1924 + (t1510 * t1541 + t1513 * t1540 + t1516 * t1539) * t1525) * pkin(1), -t1542 * t1508 - t1544 * t1505 - t1546 * t1502 + (-t1392 * t1644 - t1393 * t1640 - t1394 * t1620 + (t1392 * t1598 + t1393 * t1597 + t1394 * t1596) * t1774) * t1525, -t1543 * t1508 - t1545 * t1505 - t1547 * t1502 + (-t1400 * t1620 / 0.2e1 - t1399 * t1640 / 0.2e1 - t1398 * t1644 / 0.2e1 + (t1398 * t1598 + t1399 * t1597 + t1400 * t1596) * t1522) * t1525, t1508 * t1567 + t1505 * t1569 + t1502 * t1571 + (t1415 * t1617 + t1414 * t1618 + t1413 * t1619 + (t1413 * t1658 + t1414 * t1656 + t1415 * t1654) * t1522) * t1525, -t1508 * t1566 - t1505 * t1568 - t1502 * t1570 + (t1304 * t1412 * t1663 + t1303 * t1411 * t1664 + t1302 * t1410 * t1665 + (-t1410 * t1658 - t1411 * t1656 - t1412 * t1654) * t1522) * t1525, (t1277 * t1663 + t1276 * t1664 + t1275 * t1665 + (-t1275 * t1707 - t1276 * t1706 - t1277 * t1705) * t1522) * t1525, t1515 * t1280 * t1663 + t1512 * t1279 * t1664 + t1509 * t1278 * t1665 + t1508 * t1527 + t1505 * t1529 + t1502 * t1528 + (t1314 * t1665 + t1316 * t1664 + t1318 * t1663) * t1525 + (-t1679 * t1830 - t1677 * t1828 - t1675 * t1826 + (-t1332 * t1707 - t1333 * t1706 - t1334 * t1705) * t1525) * t1522, -t1424 * t1280 / 0.2e1 - t1423 * t1279 / 0.2e1 - t1422 * t1278 / 0.2e1 + t1508 * t1530 + t1505 * t1531 + t1502 * t1532 + (t1315 * t1665 + t1317 * t1664 + t1319 * t1663) * t1525 + (t1679 + t1677 + t1675 + (-t1329 * t1707 - t1330 * t1706 - t1331 * t1705) * t1525) * t1522, 0;];
tau_reg  = t1;