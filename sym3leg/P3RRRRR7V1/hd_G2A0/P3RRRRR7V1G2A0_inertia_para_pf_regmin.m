% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRRRR7V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x18]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:52
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:45:38
% EndTime: 2020-08-07 03:45:55
% DurationCPUTime: 19.81s
% Computational Cost: add. (18447->706), mult. (32532->1317), div. (3480->27), fcn. (26799->60), ass. (0->589)
t1829 = -2 * pkin(1);
t1828 = 2 * pkin(1);
t1827 = 2 * pkin(2);
t1485 = (pkin(4) + pkin(5));
t1826 = 2 * t1485;
t1496 = 1 / pkin(1);
t1825 = 2 * t1496;
t1824 = pkin(4) / 0.2e1;
t1468 = sin(qJ(2,3));
t1823 = pkin(4) * t1468;
t1471 = sin(qJ(2,2));
t1822 = pkin(4) * t1471;
t1474 = sin(qJ(2,1));
t1821 = pkin(4) * t1474;
t1476 = cos(qJ(3,3));
t1820 = t1476 * pkin(2);
t1479 = cos(qJ(3,2));
t1819 = t1479 * pkin(2);
t1482 = cos(qJ(3,1));
t1818 = t1482 * pkin(2);
t1817 = -qJ(3,1) + qJ(1,1);
t1816 = qJ(3,1) + qJ(1,1);
t1815 = -qJ(3,2) + qJ(1,2);
t1814 = qJ(3,2) + qJ(1,2);
t1813 = -qJ(3,3) + qJ(1,3);
t1812 = qJ(3,3) + qJ(1,3);
t1811 = 0.2e1 * qJ(2,1) + qJ(1,1);
t1810 = -0.2e1 * qJ(2,1) + qJ(1,1);
t1809 = 0.2e1 * qJ(2,2) + qJ(1,2);
t1808 = -0.2e1 * qJ(2,2) + qJ(1,2);
t1807 = 0.2e1 * qJ(2,3) + qJ(1,3);
t1806 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1478 = cos(qJ(1,3));
t1486 = 0.2e1 * qJ(3,3);
t1345 = (sin(qJ(2,3) - t1813) - sin(qJ(2,3) + t1812)) * t1826 + (-cos(0.2e1 * qJ(3,3) - t1806) - cos(t1486 + t1807) - 0.2e1 * t1478) * pkin(2) + (-cos(qJ(3,3) - t1806) - cos(qJ(3,3) + t1807) - cos(t1813) - cos(t1812)) * pkin(1);
t1461 = qJ(2,3) + qJ(3,3);
t1374 = (-sin(t1486 + qJ(2,3)) + t1468) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t1461)) * pkin(1);
t1805 = t1345 / t1374;
t1481 = cos(qJ(1,2));
t1489 = 0.2e1 * qJ(3,2);
t1346 = (sin(qJ(2,2) - t1815) - sin(qJ(2,2) + t1814)) * t1826 + (-cos(0.2e1 * qJ(3,2) - t1808) - cos(t1489 + t1809) - 0.2e1 * t1481) * pkin(2) + (-cos(qJ(3,2) - t1808) - cos(qJ(3,2) + t1809) - cos(t1815) - cos(t1814)) * pkin(1);
t1462 = qJ(2,2) + qJ(3,2);
t1375 = (-sin(t1489 + qJ(2,2)) + t1471) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t1462)) * pkin(1);
t1804 = t1346 / t1375;
t1484 = cos(qJ(1,1));
t1492 = 0.2e1 * qJ(3,1);
t1347 = (sin(qJ(2,1) - t1817) - sin(qJ(2,1) + t1816)) * t1826 + (-cos(0.2e1 * qJ(3,1) - t1810) - cos(t1492 + t1811) - 0.2e1 * t1484) * pkin(2) + (-cos(qJ(3,1) - t1810) - cos(qJ(3,1) + t1811) - cos(t1817) - cos(t1816)) * pkin(1);
t1463 = qJ(2,1) + qJ(3,1);
t1376 = (-sin(t1492 + qJ(2,1)) + t1474) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t1463)) * pkin(1);
t1803 = t1347 / t1376;
t1469 = sin(qJ(1,3));
t1423 = pkin(1) + t1820;
t1477 = cos(qJ(2,3));
t1467 = sin(qJ(3,3));
t1757 = t1467 * t1468;
t1545 = pkin(2) * t1757 - t1423 * t1477;
t1743 = t1478 * t1485;
t1356 = -t1469 * t1545 - t1743;
t1756 = t1467 * t1477;
t1392 = pkin(2) * t1756 + t1468 * t1423;
t1464 = legFrame(3,2);
t1434 = sin(t1464);
t1437 = cos(t1464);
t1350 = t1356 * t1434 - t1392 * t1437;
t1440 = 0.1e1 / t1467;
t1802 = t1350 * t1440;
t1351 = -t1356 * t1437 - t1392 * t1434;
t1801 = t1351 * t1440;
t1472 = sin(qJ(1,2));
t1425 = pkin(1) + t1819;
t1480 = cos(qJ(2,2));
t1470 = sin(qJ(3,2));
t1753 = t1470 * t1471;
t1544 = pkin(2) * t1753 - t1425 * t1480;
t1740 = t1481 * t1485;
t1357 = -t1472 * t1544 - t1740;
t1752 = t1470 * t1480;
t1393 = pkin(2) * t1752 + t1471 * t1425;
t1465 = legFrame(2,2);
t1435 = sin(t1465);
t1438 = cos(t1465);
t1352 = t1357 * t1435 - t1393 * t1438;
t1444 = 0.1e1 / t1470;
t1800 = t1352 * t1444;
t1353 = -t1357 * t1438 - t1393 * t1435;
t1799 = t1353 * t1444;
t1475 = sin(qJ(1,1));
t1427 = pkin(1) + t1818;
t1483 = cos(qJ(2,1));
t1473 = sin(qJ(3,1));
t1749 = t1473 * t1474;
t1543 = pkin(2) * t1749 - t1427 * t1483;
t1737 = t1484 * t1485;
t1358 = -t1475 * t1543 - t1737;
t1748 = t1473 * t1483;
t1394 = pkin(2) * t1748 + t1474 * t1427;
t1466 = legFrame(1,2);
t1436 = sin(t1466);
t1439 = cos(t1466);
t1354 = t1358 * t1436 - t1394 * t1439;
t1448 = 0.1e1 / t1473;
t1798 = t1354 * t1448;
t1355 = -t1358 * t1439 - t1394 * t1436;
t1797 = t1355 * t1448;
t1362 = -t1469 * t1485 + t1478 * t1545;
t1796 = t1362 * t1440;
t1363 = -t1472 * t1485 + t1481 * t1544;
t1795 = t1363 * t1444;
t1364 = -t1475 * t1485 + t1484 * t1543;
t1794 = t1364 * t1448;
t1386 = 0.1e1 / t1545;
t1793 = t1386 * t1440;
t1792 = 0.1e1 / t1545 ^ 2 / t1467 ^ 2;
t1388 = 0.1e1 / t1544;
t1791 = t1388 * t1444;
t1790 = 0.1e1 / t1544 ^ 2 / t1470 ^ 2;
t1390 = 0.1e1 / t1543;
t1789 = t1390 * t1448;
t1788 = 0.1e1 / t1543 ^ 2 / t1473 ^ 2;
t1416 = pkin(2) * cos(t1461) + t1477 * pkin(1);
t1410 = 0.1e1 / t1416;
t1787 = t1410 * t1469;
t1786 = t1410 * t1478;
t1411 = 0.1e1 / t1416 ^ 2;
t1454 = t1478 ^ 2;
t1785 = t1411 * t1454;
t1417 = pkin(2) * cos(t1462) + t1480 * pkin(1);
t1412 = 0.1e1 / t1417;
t1784 = t1412 * t1472;
t1783 = t1412 * t1481;
t1413 = 0.1e1 / t1417 ^ 2;
t1457 = t1481 ^ 2;
t1782 = t1413 * t1457;
t1418 = pkin(2) * cos(t1463) + t1483 * pkin(1);
t1414 = 0.1e1 / t1418;
t1781 = t1414 * t1475;
t1780 = t1414 * t1484;
t1415 = 0.1e1 / t1418 ^ 2;
t1460 = t1484 ^ 2;
t1779 = t1415 * t1460;
t1452 = t1476 ^ 2;
t1419 = pkin(1) * t1476 + t1452 * t1827 - pkin(2);
t1778 = t1419 * t1468;
t1777 = t1419 * t1469;
t1455 = t1479 ^ 2;
t1420 = pkin(1) * t1479 + t1455 * t1827 - pkin(2);
t1776 = t1420 * t1471;
t1775 = t1420 * t1472;
t1458 = t1482 ^ 2;
t1421 = t1482 * pkin(1) + t1458 * t1827 - pkin(2);
t1774 = t1421 * t1474;
t1773 = t1421 * t1475;
t1772 = (pkin(1) + 0.2e1 * t1820) * t1467;
t1771 = (pkin(1) + 0.2e1 * t1819) * t1470;
t1770 = (pkin(1) + 0.2e1 * t1818) * t1473;
t1769 = t1434 * t1478;
t1768 = t1435 * t1481;
t1767 = t1436 * t1484;
t1766 = t1437 * t1478;
t1765 = t1438 * t1481;
t1764 = t1439 * t1484;
t1763 = t1440 * t1469;
t1443 = t1469 ^ 2;
t1762 = t1443 * t1411;
t1761 = t1444 * t1472;
t1447 = t1472 ^ 2;
t1760 = t1447 * t1413;
t1759 = t1448 * t1475;
t1451 = t1475 ^ 2;
t1758 = t1451 * t1415;
t1755 = t1468 * t1469;
t1754 = t1468 * t1477;
t1751 = t1471 * t1472;
t1750 = t1471 * t1480;
t1747 = t1474 * t1475;
t1746 = t1474 * t1483;
t1745 = t1476 * t1467;
t1744 = t1476 * t1477;
t1742 = t1479 * t1470;
t1741 = t1479 * t1480;
t1739 = t1482 * t1473;
t1738 = t1482 * t1483;
t1495 = 1 / pkin(2);
t1736 = t1495 * t1496;
t1735 = -0.2e1 * t1756;
t1734 = -0.2e1 * t1752;
t1733 = -0.2e1 * t1748;
t1732 = -0.2e1 * t1744;
t1731 = -0.2e1 * t1741;
t1730 = -0.2e1 * t1738;
t1729 = pkin(2) * t1745;
t1728 = pkin(2) * t1742;
t1727 = pkin(2) * t1739;
t1726 = pkin(4) * t1410 * t1477;
t1725 = pkin(4) * t1412 * t1480;
t1724 = pkin(4) * t1414 * t1483;
t1377 = -t1743 * t1757 + (t1452 - 0.1e1) * t1469 * pkin(2);
t1453 = t1477 ^ 2;
t1671 = t1467 * t1755;
t1518 = pkin(1) * t1671 + (t1671 * t1827 + t1743) * t1476;
t1321 = (t1434 * t1772 + t1437 * t1777) * t1453 + (t1434 * t1778 - t1437 * t1518) * t1477 - t1377 * t1437 - t1434 * t1729;
t1711 = t1496 * t1793;
t1301 = t1321 * t1711;
t1674 = t1440 * t1736;
t1270 = t1351 * t1674 - t1301;
t1723 = t1270 * t1787;
t1378 = -t1740 * t1753 + (t1455 - 0.1e1) * t1472 * pkin(2);
t1456 = t1480 ^ 2;
t1668 = t1470 * t1751;
t1517 = pkin(1) * t1668 + (t1668 * t1827 + t1740) * t1479;
t1323 = (t1435 * t1771 + t1438 * t1775) * t1456 + (t1435 * t1776 - t1438 * t1517) * t1480 - t1378 * t1438 - t1435 * t1728;
t1709 = t1496 * t1791;
t1303 = t1323 * t1709;
t1673 = t1444 * t1736;
t1272 = t1353 * t1673 - t1303;
t1722 = t1272 * t1784;
t1320 = (-t1434 * t1777 + t1437 * t1772) * t1453 + (t1434 * t1518 + t1437 * t1778) * t1477 + t1377 * t1434 - t1437 * t1729;
t1721 = t1320 * t1793;
t1720 = t1321 * t1793;
t1322 = (-t1435 * t1775 + t1438 * t1771) * t1456 + (t1435 * t1517 + t1438 * t1776) * t1480 + t1378 * t1435 - t1438 * t1728;
t1719 = t1322 * t1791;
t1718 = t1323 * t1791;
t1379 = -t1737 * t1749 + (t1458 - 0.1e1) * t1475 * pkin(2);
t1459 = t1483 ^ 2;
t1665 = t1473 * t1747;
t1516 = pkin(1) * t1665 + (t1665 * t1827 + t1737) * t1482;
t1324 = (-t1436 * t1773 + t1439 * t1770) * t1459 + (t1436 * t1516 + t1439 * t1774) * t1483 + t1379 * t1436 - t1439 * t1727;
t1717 = t1324 * t1789;
t1325 = (t1436 * t1770 + t1439 * t1773) * t1459 + (t1436 * t1774 - t1439 * t1516) * t1483 - t1379 * t1439 - t1436 * t1727;
t1716 = t1325 * t1789;
t1663 = t1805 / 0.2e1;
t1335 = t1496 * t1663;
t1326 = t1362 * t1674 + t1335;
t1715 = t1326 * t1793;
t1662 = t1804 / 0.2e1;
t1336 = t1496 * t1662;
t1327 = t1363 * t1673 + t1336;
t1714 = t1327 * t1791;
t1661 = t1803 / 0.2e1;
t1337 = t1496 * t1661;
t1672 = t1448 * t1736;
t1328 = t1364 * t1672 + t1337;
t1713 = t1328 * t1789;
t1712 = t1386 * t1763;
t1710 = t1388 * t1761;
t1708 = t1390 * t1759;
t1707 = t1496 * t1789;
t1406 = t1738 - t1749;
t1706 = t1406 * t1781;
t1409 = t1482 * t1474 + t1748;
t1705 = t1409 * t1781;
t1704 = t1410 * t1769;
t1703 = t1410 * t1766;
t1702 = t1410 * t1763;
t1701 = t1440 * t1786;
t1700 = t1410 * t1755;
t1699 = t1468 * t1786;
t1698 = t1411 * t1754;
t1697 = t1411 * t1469 * t1478;
t1696 = t1412 * t1768;
t1695 = t1412 * t1765;
t1694 = t1412 * t1761;
t1693 = t1444 * t1783;
t1692 = t1412 * t1751;
t1691 = t1471 * t1783;
t1690 = t1413 * t1750;
t1689 = t1413 * t1472 * t1481;
t1688 = t1414 * t1767;
t1687 = t1414 * t1764;
t1686 = t1414 * t1759;
t1685 = t1448 * t1780;
t1684 = t1414 * t1747;
t1683 = t1474 * t1780;
t1682 = t1415 * t1746;
t1681 = t1415 * t1475 * t1484;
t1428 = t1434 ^ 2;
t1680 = t1428 * t1785;
t1429 = t1435 ^ 2;
t1679 = t1429 * t1782;
t1430 = t1436 ^ 2;
t1678 = t1430 * t1779;
t1431 = t1437 ^ 2;
t1677 = t1431 * t1785;
t1432 = t1438 ^ 2;
t1676 = t1432 * t1782;
t1433 = t1439 ^ 2;
t1675 = t1433 * t1779;
t1300 = t1320 * t1711;
t1269 = t1350 * t1674 - t1300;
t1670 = t1269 * t1787;
t1669 = t1326 * t1787;
t1302 = t1322 * t1709;
t1271 = t1352 * t1673 - t1302;
t1667 = t1271 * t1784;
t1666 = t1327 * t1784;
t1664 = t1328 * t1781;
t1660 = t1736 / 0.2e1;
t1659 = t1410 * t1453 * t1828;
t1658 = t1412 * t1456 * t1828;
t1657 = t1414 * t1459 * t1828;
t1656 = t1469 * t1726;
t1655 = t1478 * t1726;
t1654 = t1472 * t1725;
t1653 = t1481 * t1725;
t1652 = t1475 * t1724;
t1651 = t1484 * t1724;
t1650 = pkin(4) * t1700;
t1649 = pkin(4) * t1692;
t1648 = pkin(4) * t1684;
t1647 = t1787 * t1805;
t1646 = t1784 * t1804;
t1645 = t1781 * t1803;
t1404 = t1744 - t1757;
t1644 = t1404 * t1704;
t1643 = t1404 * t1703;
t1642 = t1404 * t1702;
t1405 = t1741 - t1753;
t1641 = t1405 * t1696;
t1640 = t1405 * t1695;
t1639 = t1405 * t1694;
t1638 = t1406 * t1688;
t1637 = t1406 * t1687;
t1636 = t1406 * t1686;
t1407 = t1476 * t1468 + t1756;
t1635 = t1407 * t1704;
t1634 = t1407 * t1703;
t1633 = t1407 * t1702;
t1408 = t1479 * t1471 + t1752;
t1632 = t1408 * t1696;
t1631 = t1408 * t1695;
t1630 = t1408 * t1694;
t1629 = t1409 * t1688;
t1628 = t1409 * t1687;
t1627 = t1409 * t1686;
t1626 = t1434 * t1701;
t1625 = t1434 * t1699;
t1624 = t1437 * t1701;
t1623 = t1437 * t1699;
t1622 = t1454 * t1698;
t1621 = t1435 * t1693;
t1620 = t1435 * t1691;
t1619 = t1438 * t1693;
t1618 = t1438 * t1691;
t1617 = t1457 * t1690;
t1616 = t1436 * t1685;
t1615 = t1436 * t1683;
t1614 = t1439 * t1685;
t1613 = t1439 * t1683;
t1612 = t1460 * t1682;
t1611 = t1434 * t1697;
t1610 = t1435 * t1689;
t1609 = t1436 * t1681;
t1608 = t1437 * t1434 * t1785;
t1607 = t1437 * t1697;
t1606 = t1438 * t1435 * t1782;
t1605 = t1438 * t1689;
t1604 = t1439 * t1436 * t1779;
t1603 = t1439 * t1681;
t1602 = t1440 * t1660;
t1601 = t1444 * t1660;
t1600 = t1448 * t1660;
t1599 = t1478 * t1659;
t1598 = t1481 * t1658;
t1597 = t1484 * t1657;
t1596 = t1467 * t1655;
t1595 = t1476 * t1655;
t1594 = t1470 * t1653;
t1593 = t1479 * t1653;
t1592 = t1473 * t1651;
t1591 = t1482 * t1651;
t1590 = pkin(4) * t1625;
t1589 = pkin(4) * t1623;
t1588 = pkin(4) * t1620;
t1587 = pkin(4) * t1618;
t1586 = pkin(4) * t1615;
t1585 = pkin(4) * t1613;
t1584 = t1386 * t1642;
t1583 = t1386 * t1633;
t1582 = t1388 * t1639;
t1581 = t1388 * t1630;
t1580 = t1390 * t1636;
t1579 = t1390 * t1627;
t1578 = t1404 * t1626;
t1577 = t1404 * t1624;
t1576 = t1405 * t1621;
t1575 = t1405 * t1619;
t1574 = t1406 * t1616;
t1573 = t1406 * t1614;
t1572 = t1407 * t1626;
t1571 = t1407 * t1624;
t1570 = t1408 * t1621;
t1569 = t1408 * t1619;
t1568 = t1409 * t1616;
t1567 = t1409 * t1614;
t1566 = t1697 * t1754;
t1565 = t1689 * t1750;
t1564 = t1681 * t1746;
t1563 = t1663 * t1793;
t1562 = -t1647 / 0.2e1;
t1561 = t1662 * t1791;
t1560 = -t1646 / 0.2e1;
t1559 = t1661 * t1789;
t1558 = -t1645 / 0.2e1;
t1557 = -t1769 * t1805 / 0.2e1;
t1556 = -t1768 * t1804 / 0.2e1;
t1555 = -t1767 * t1803 / 0.2e1;
t1554 = t1663 * t1766;
t1553 = t1662 * t1765;
t1552 = t1661 * t1764;
t1551 = t1437 * t1596;
t1550 = t1438 * t1594;
t1549 = t1439 * t1592;
t1548 = t1437 * t1595;
t1547 = t1438 * t1593;
t1546 = t1439 * t1591;
t1542 = t1320 * t1386 * t1626;
t1541 = t1321 * t1386 * t1624;
t1540 = t1322 * t1388 * t1621;
t1539 = t1323 * t1388 * t1619;
t1538 = t1324 * t1390 * t1616;
t1537 = t1325 * t1390 * t1614;
t1536 = t1386 * t1578;
t1535 = t1386 * t1577;
t1534 = t1386 * t1572;
t1533 = t1386 * t1571;
t1532 = t1388 * t1576;
t1531 = t1388 * t1575;
t1530 = t1388 * t1570;
t1529 = t1388 * t1569;
t1528 = t1390 * t1574;
t1527 = t1390 * t1573;
t1526 = t1390 * t1568;
t1525 = t1390 * t1567;
t1524 = t1410 * t1557;
t1523 = t1414 * t1555;
t1522 = t1412 * t1556;
t1521 = t1410 * t1554;
t1520 = t1412 * t1553;
t1519 = t1414 * t1552;
t1515 = t1326 * t1823 + t1469 * t1659;
t1514 = t1327 * t1822 + t1472 * t1658;
t1513 = t1328 * t1821 + t1475 * t1657;
t1512 = t1269 * t1823 + t1434 * t1599;
t1511 = -t1270 * t1823 + t1437 * t1599;
t1510 = t1271 * t1822 + t1435 * t1598;
t1509 = -t1272 * t1822 + t1438 * t1598;
t1304 = t1324 * t1707;
t1273 = t1354 * t1672 - t1304;
t1508 = t1273 * t1821 + t1436 * t1597;
t1305 = t1325 * t1707;
t1274 = t1355 * t1672 - t1305;
t1507 = -t1274 * t1821 + t1439 * t1597;
t1506 = t1386 * (t1320 * t1437 - t1321 * t1434) * t1701;
t1505 = t1388 * (t1322 * t1438 - t1323 * t1435) * t1693;
t1504 = t1390 * (t1324 * t1439 - t1325 * t1436) * t1685;
t1503 = t1410 * (t1320 * t1712 + t1557);
t1502 = t1410 * (t1321 * t1712 + t1554);
t1501 = t1412 * (t1322 * t1710 + t1556);
t1500 = t1412 * (t1323 * t1710 + t1553);
t1499 = t1414 * (t1324 * t1708 + t1555);
t1498 = t1414 * (t1325 * t1708 + t1552);
t1497 = 1 / pkin(1) ^ 2;
t1450 = t1474 ^ 2;
t1446 = t1471 ^ 2;
t1442 = t1468 ^ 2;
t1403 = t1409 ^ 2;
t1402 = t1408 ^ 2;
t1401 = t1407 ^ 2;
t1385 = t1482 * t1652;
t1384 = t1479 * t1654;
t1383 = t1476 * t1656;
t1382 = t1473 * t1652;
t1381 = t1470 * t1654;
t1380 = t1467 * t1656;
t1370 = t1436 * t1591;
t1369 = t1435 * t1593;
t1368 = t1434 * t1595;
t1367 = t1436 * t1592;
t1366 = t1435 * t1594;
t1365 = t1434 * t1596;
t1361 = (t1458 - 0.1e1 / 0.2e1) * t1746 + (t1459 - 0.1e1 / 0.2e1) * t1739;
t1360 = (t1455 - 0.1e1 / 0.2e1) * t1750 + (t1456 - 0.1e1 / 0.2e1) * t1742;
t1359 = (t1452 - 0.1e1 / 0.2e1) * t1754 + (t1453 - 0.1e1 / 0.2e1) * t1745;
t1349 = -t1603 - t1605 - t1607;
t1348 = t1609 + t1610 + t1611;
t1344 = -t1604 - t1606 - t1608;
t1343 = -t1442 * t1607 - t1446 * t1605 - t1450 * t1603;
t1342 = t1442 * t1611 + t1446 * t1610 + t1450 * t1609;
t1341 = -t1442 * t1608 - t1446 * t1606 - t1450 * t1604;
t1340 = -0.2e1 * t1437 * t1566 - 0.2e1 * t1438 * t1565 - 0.2e1 * t1439 * t1564;
t1339 = 0.2e1 * t1434 * t1566 + 0.2e1 * t1435 * t1565 + 0.2e1 * t1436 * t1564;
t1338 = -0.2e1 * t1604 * t1746 - 0.2e1 * t1606 * t1750 - 0.2e1 * t1608 * t1754;
t1334 = -t1401 * t1607 - t1402 * t1605 - t1403 * t1603;
t1333 = t1401 * t1611 + t1402 * t1610 + t1403 * t1609;
t1332 = -t1401 * t1608 - t1402 * t1606 - t1403 * t1604;
t1331 = t1648 + t1661;
t1330 = t1649 + t1662;
t1329 = t1650 + t1663;
t1319 = -t1473 * t1331 + t1385;
t1318 = -t1470 * t1330 + t1384;
t1317 = -t1467 * t1329 + t1383;
t1316 = t1482 * t1331 + t1382;
t1315 = t1479 * t1330 + t1381;
t1314 = t1476 * t1329 + t1380;
t1313 = -pkin(1) * t1684 + t1328 * t1824;
t1312 = t1648 + (t1364 * t1600 + t1337) * t1828;
t1311 = -pkin(1) * t1692 + t1327 * t1824;
t1310 = t1649 + (t1363 * t1601 + t1336) * t1828;
t1309 = -pkin(1) * t1700 + t1326 * t1824;
t1308 = t1650 + (t1362 * t1602 + t1335) * t1828;
t1307 = -0.4e1 * t1359 * t1607 - 0.4e1 * t1360 * t1605 - 0.4e1 * t1361 * t1603;
t1306 = 0.4e1 * t1359 * t1611 + 0.4e1 * t1360 * t1610 + 0.4e1 * t1361 * t1609;
t1299 = -0.4e1 * t1359 * t1608 - 0.4e1 * t1360 * t1606 - 0.4e1 * t1361 * t1604;
t1298 = -t1312 * t1473 + t1385;
t1297 = t1312 * t1482 + t1382;
t1296 = -t1310 * t1470 + t1384;
t1295 = t1310 * t1479 + t1381;
t1294 = -t1308 * t1467 + t1383;
t1293 = t1308 * t1476 + t1380;
t1292 = -t1585 - t1716;
t1291 = t1586 - t1717;
t1290 = -t1587 - t1718;
t1289 = t1588 - t1719;
t1288 = -t1589 - t1720;
t1287 = t1590 - t1721;
t1286 = -t1473 * t1292 - t1546;
t1285 = -t1473 * t1291 + t1370;
t1284 = -t1470 * t1290 - t1547;
t1283 = -t1470 * t1289 + t1369;
t1282 = -t1467 * t1288 - t1548;
t1281 = -t1467 * t1287 + t1368;
t1280 = t1482 * t1292 - t1549;
t1279 = t1482 * t1291 + t1367;
t1278 = t1479 * t1290 - t1550;
t1277 = t1479 * t1289 + t1366;
t1276 = t1476 * t1288 - t1551;
t1275 = t1476 * t1287 + t1365;
t1268 = pkin(1) * t1613 + t1274 * t1824;
t1267 = -pkin(1) * t1615 + t1273 * t1824;
t1266 = t1585 + (t1355 * t1600 - t1305) * t1829;
t1265 = t1586 + (t1354 * t1600 - t1304) * t1828;
t1264 = pkin(1) * t1618 + t1272 * t1824;
t1263 = -pkin(1) * t1620 + t1271 * t1824;
t1262 = t1587 + (t1353 * t1601 - t1303) * t1829;
t1261 = t1588 + (t1352 * t1601 - t1302) * t1828;
t1260 = pkin(1) * t1623 + t1270 * t1824;
t1259 = -pkin(1) * t1625 + t1269 * t1824;
t1258 = t1589 + (t1351 * t1602 - t1301) * t1829;
t1257 = t1590 + (t1350 * t1602 - t1300) * t1828;
t1256 = -t1266 * t1482 - t1549;
t1255 = t1266 * t1473 - t1546;
t1254 = -t1265 * t1473 + t1370;
t1253 = t1265 * t1482 + t1367;
t1252 = -t1262 * t1479 - t1550;
t1251 = t1262 * t1470 - t1547;
t1250 = -t1261 * t1470 + t1369;
t1249 = t1261 * t1479 + t1366;
t1248 = -t1258 * t1476 - t1551;
t1247 = t1258 * t1467 - t1548;
t1246 = -t1257 * t1467 + t1368;
t1245 = t1257 * t1476 + t1365;
t1244 = t1313 * t1730 + t1473 * t1513;
t1243 = t1313 * t1733 - t1482 * t1513;
t1242 = t1311 * t1731 + t1470 * t1514;
t1241 = t1311 * t1734 - t1479 * t1514;
t1240 = t1309 * t1732 + t1467 * t1515;
t1239 = t1309 * t1735 - t1476 * t1515;
t1238 = t1268 * t1730 - t1473 * t1507;
t1237 = t1267 * t1730 + t1473 * t1508;
t1236 = t1268 * t1733 + t1482 * t1507;
t1235 = t1267 * t1733 - t1482 * t1508;
t1234 = t1264 * t1731 - t1470 * t1509;
t1233 = t1263 * t1731 + t1470 * t1510;
t1232 = t1264 * t1734 + t1479 * t1509;
t1231 = t1263 * t1734 - t1479 * t1510;
t1230 = t1260 * t1732 - t1467 * t1511;
t1229 = t1259 * t1732 + t1467 * t1512;
t1228 = t1260 * t1735 + t1476 * t1511;
t1227 = t1259 * t1735 - t1476 * t1512;
t1226 = (-t1321 * t1563 - t1323 * t1561 - t1325 * t1559) * t1497;
t1225 = (-t1320 * t1563 - t1322 * t1561 - t1324 * t1559) * t1497;
t1224 = (t1320 * t1321 * t1792 + t1322 * t1323 * t1790 + t1324 * t1325 * t1788) * t1497;
t1223 = (t1477 * t1502 + t1480 * t1500 + t1483 * t1498) * t1496;
t1222 = (t1468 * t1502 + t1471 * t1500 + t1474 * t1498) * t1496;
t1221 = (t1477 * t1503 + t1480 * t1501 + t1483 * t1499) * t1496;
t1220 = (t1468 * t1503 + t1471 * t1501 + t1474 * t1499) * t1496;
t1219 = (-t1477 * t1506 - t1480 * t1505 - t1483 * t1504) * t1496;
t1218 = (-t1468 * t1506 - t1471 * t1505 - t1474 * t1504) * t1496;
t1 = [t1675 + t1676 + t1677, 0, 0, t1442 * t1677 + t1446 * t1676 + t1450 * t1675, 0.2e1 * t1431 * t1622 + 0.2e1 * t1432 * t1617 + 0.2e1 * t1433 * t1612, (-t1468 * t1541 - t1471 * t1539 - t1474 * t1537) * t1825, (-t1477 * t1541 - t1480 * t1539 - t1483 * t1537) * t1825, (t1321 ^ 2 * t1792 + t1323 ^ 2 * t1790 + t1325 ^ 2 * t1788) * t1497, 0, 0, t1401 * t1677 + t1402 * t1676 + t1403 * t1675, 0.4e1 * t1359 * t1677 + 0.4e1 * t1360 * t1676 + 0.4e1 * t1361 * t1675, t1270 * t1634 + t1272 * t1631 + t1274 * t1628 + (-t1321 * t1533 - t1323 * t1529 - t1325 * t1525 + (t1351 * t1571 + t1353 * t1569 + t1355 * t1567) * t1495) * t1496, t1270 * t1643 + t1272 * t1640 + t1274 * t1637 + (-t1321 * t1535 - t1323 * t1531 - t1325 * t1527 + (t1351 * t1577 + t1353 * t1575 + t1355 * t1573) * t1495) * t1496, (-t1270 * t1720 - t1272 * t1718 - t1274 * t1716 + (t1270 * t1801 + t1272 * t1799 + t1274 * t1797) * t1495) * t1496, t1228 * t1703 + t1232 * t1695 + t1236 * t1687 + (-t1248 * t1720 - t1252 * t1718 - t1256 * t1716 + (t1276 * t1801 + t1278 * t1799 + t1280 * t1797) * t1495) * t1496, t1230 * t1703 + t1234 * t1695 + t1238 * t1687 + (-t1247 * t1720 - t1251 * t1718 - t1255 * t1716 + (t1282 * t1801 + t1284 * t1799 + t1286 * t1797) * t1495) * t1496, 1; t1344, 0, 0, t1341, t1338, t1218, t1219, t1224, 0, 0, t1332, t1299, t1269 * t1634 + t1271 * t1631 + t1273 * t1628 + (t1321 * t1534 + t1323 * t1530 + t1325 * t1526 + (-t1351 * t1572 - t1353 * t1570 - t1355 * t1568) * t1495) * t1496, t1269 * t1643 + t1271 * t1640 + t1273 * t1637 + (t1321 * t1536 + t1323 * t1532 + t1325 * t1528 + (-t1351 * t1578 - t1353 * t1576 - t1355 * t1574) * t1495) * t1496, (-t1269 * t1720 - t1271 * t1718 - t1273 * t1716 + (t1269 * t1801 + t1271 * t1799 + t1273 * t1797) * t1495) * t1496, t1227 * t1703 + t1231 * t1695 + t1235 * t1687 + (-t1245 * t1720 - t1249 * t1718 - t1253 * t1716 + (t1275 * t1801 + t1277 * t1799 + t1279 * t1797) * t1495) * t1496, t1229 * t1703 + t1233 * t1695 + t1237 * t1687 + (-t1246 * t1720 - t1250 * t1718 - t1254 * t1716 + (t1281 * t1801 + t1283 * t1799 + t1285 * t1797) * t1495) * t1496, 0; t1349, 0, 0, t1343, t1340, t1222, t1223, t1226, 0, 0, t1334, t1307, t1326 * t1634 + t1327 * t1631 + t1328 * t1628 + (t1321 * t1583 + t1323 * t1581 + t1325 * t1579 + (-t1351 * t1633 - t1353 * t1630 - t1355 * t1627) * t1495) * t1496, t1326 * t1643 + t1327 * t1640 + t1328 * t1637 + (t1321 * t1584 + t1323 * t1582 + t1325 * t1580 + (-t1351 * t1642 - t1353 * t1639 - t1355 * t1636) * t1495) * t1496, (-t1321 * t1715 - t1323 * t1714 - t1325 * t1713 + (t1326 * t1801 + t1327 * t1799 + t1328 * t1797) * t1495) * t1496, t1239 * t1703 + t1241 * t1695 + t1243 * t1687 + (-t1293 * t1720 - t1295 * t1718 - t1297 * t1716 + (t1314 * t1801 + t1315 * t1799 + t1316 * t1797) * t1495) * t1496, t1240 * t1703 + t1242 * t1695 + t1244 * t1687 + (-t1294 * t1720 - t1296 * t1718 - t1298 * t1716 + (t1317 * t1801 + t1318 * t1799 + t1319 * t1797) * t1495) * t1496, 0; t1344, 0, 0, t1341, t1338, t1218, t1219, t1224, 0, 0, t1332, t1299, -t1270 * t1635 - t1272 * t1632 - t1274 * t1629 + (-t1320 * t1533 - t1322 * t1529 - t1324 * t1525 + (t1350 * t1571 + t1352 * t1569 + t1354 * t1567) * t1495) * t1496, -t1270 * t1644 - t1272 * t1641 - t1274 * t1638 + (-t1320 * t1535 - t1322 * t1531 - t1324 * t1527 + (t1350 * t1577 + t1352 * t1575 + t1354 * t1573) * t1495) * t1496, (-t1270 * t1721 - t1272 * t1719 - t1274 * t1717 + (t1270 * t1802 + t1272 * t1800 + t1274 * t1798) * t1495) * t1496, -t1228 * t1704 - t1232 * t1696 - t1236 * t1688 + (-t1248 * t1721 - t1252 * t1719 - t1256 * t1717 + (t1276 * t1802 + t1278 * t1800 + t1280 * t1798) * t1495) * t1496, -t1230 * t1704 - t1234 * t1696 - t1238 * t1688 + (-t1247 * t1721 - t1251 * t1719 - t1255 * t1717 + (t1282 * t1802 + t1284 * t1800 + t1286 * t1798) * t1495) * t1496, 0; t1678 + t1679 + t1680, 0, 0, t1442 * t1680 + t1446 * t1679 + t1450 * t1678, 0.2e1 * t1428 * t1622 + 0.2e1 * t1429 * t1617 + 0.2e1 * t1430 * t1612, (t1468 * t1542 + t1471 * t1540 + t1474 * t1538) * t1825, (t1477 * t1542 + t1480 * t1540 + t1483 * t1538) * t1825, (t1320 ^ 2 * t1792 + t1322 ^ 2 * t1790 + t1324 ^ 2 * t1788) * t1497, 0, 0, t1401 * t1680 + t1402 * t1679 + t1403 * t1678, 0.4e1 * t1359 * t1680 + 0.4e1 * t1360 * t1679 + 0.4e1 * t1361 * t1678, -t1269 * t1635 - t1271 * t1632 - t1273 * t1629 + (t1320 * t1534 + t1322 * t1530 + t1324 * t1526 + (-t1350 * t1572 - t1352 * t1570 - t1354 * t1568) * t1495) * t1496, -t1269 * t1644 - t1271 * t1641 - t1273 * t1638 + (t1320 * t1536 + t1322 * t1532 + t1324 * t1528 + (-t1350 * t1578 - t1352 * t1576 - t1354 * t1574) * t1495) * t1496, (-t1269 * t1721 - t1271 * t1719 - t1273 * t1717 + (t1269 * t1802 + t1271 * t1800 + t1273 * t1798) * t1495) * t1496, -t1227 * t1704 - t1231 * t1696 - t1235 * t1688 + (-t1245 * t1721 - t1249 * t1719 - t1253 * t1717 + (t1275 * t1802 + t1277 * t1800 + t1279 * t1798) * t1495) * t1496, -t1229 * t1704 - t1233 * t1696 - t1237 * t1688 + (-t1246 * t1721 - t1250 * t1719 - t1254 * t1717 + (t1281 * t1802 + t1283 * t1800 + t1285 * t1798) * t1495) * t1496, 1; t1348, 0, 0, t1342, t1339, t1220, t1221, t1225, 0, 0, t1333, t1306, -t1326 * t1635 - t1327 * t1632 - t1328 * t1629 + (t1320 * t1583 + t1322 * t1581 + t1324 * t1579 + (-t1350 * t1633 - t1352 * t1630 - t1354 * t1627) * t1495) * t1496, -t1326 * t1644 - t1327 * t1641 - t1328 * t1638 + (t1320 * t1584 + t1322 * t1582 + t1324 * t1580 + (-t1350 * t1642 - t1352 * t1639 - t1354 * t1636) * t1495) * t1496, (-t1320 * t1715 - t1322 * t1714 - t1324 * t1713 + (t1326 * t1802 + t1327 * t1800 + t1328 * t1798) * t1495) * t1496, -t1239 * t1704 - t1241 * t1696 - t1243 * t1688 + (-t1293 * t1721 - t1295 * t1719 - t1297 * t1717 + (t1314 * t1802 + t1315 * t1800 + t1316 * t1798) * t1495) * t1496, -t1240 * t1704 - t1242 * t1696 - t1244 * t1688 + (-t1294 * t1721 - t1296 * t1719 - t1298 * t1717 + (t1317 * t1802 + t1318 * t1800 + t1319 * t1798) * t1495) * t1496, 0; t1349, 0, 0, t1343, t1340, t1222, t1223, t1226, 0, 0, t1334, t1307, -t1407 * t1723 - t1408 * t1722 - t1274 * t1705 + (t1409 * t1519 + t1408 * t1520 + t1407 * t1521 + (t1362 * t1571 + t1363 * t1569 + t1364 * t1567) * t1495) * t1496, -t1404 * t1723 - t1405 * t1722 - t1274 * t1706 + (t1406 * t1519 + t1405 * t1520 + t1404 * t1521 + (t1362 * t1577 + t1363 * t1575 + t1364 * t1573) * t1495) * t1496, (t1274 * t1661 + t1272 * t1662 + t1270 * t1663 + (t1270 * t1796 + t1272 * t1795 + t1274 * t1794) * t1495) * t1496, -t1228 * t1787 - t1232 * t1784 - t1236 * t1781 + (t1256 * t1661 + t1252 * t1662 + t1248 * t1663 + (t1276 * t1796 + t1278 * t1795 + t1280 * t1794) * t1495) * t1496, -t1230 * t1787 - t1234 * t1784 - t1238 * t1781 + (t1255 * t1661 + t1251 * t1662 + t1247 * t1663 + (t1282 * t1796 + t1284 * t1795 + t1286 * t1794) * t1495) * t1496, 0; t1348, 0, 0, t1342, t1339, t1220, t1221, t1225, 0, 0, t1333, t1306, -t1407 * t1670 - t1408 * t1667 - t1273 * t1705 + (t1409 * t1523 + t1408 * t1522 + t1407 * t1524 + (-t1362 * t1572 - t1363 * t1570 - t1364 * t1568) * t1495) * t1496, -t1404 * t1670 - t1405 * t1667 - t1273 * t1706 + (t1406 * t1523 + t1405 * t1522 + t1404 * t1524 + (-t1362 * t1578 - t1363 * t1576 - t1364 * t1574) * t1495) * t1496, (t1273 * t1661 + t1271 * t1662 + t1269 * t1663 + (t1269 * t1796 + t1271 * t1795 + t1273 * t1794) * t1495) * t1496, -t1227 * t1787 - t1231 * t1784 - t1235 * t1781 + (t1253 * t1661 + t1249 * t1662 + t1245 * t1663 + (t1275 * t1796 + t1277 * t1795 + t1279 * t1794) * t1495) * t1496, -t1229 * t1787 - t1233 * t1784 - t1237 * t1781 + (t1254 * t1661 + t1250 * t1662 + t1246 * t1663 + (t1281 * t1796 + t1283 * t1795 + t1285 * t1794) * t1495) * t1496, 0; t1758 + t1760 + t1762, 0, 0, t1442 * t1762 + t1446 * t1760 + t1450 * t1758, 0.2e1 * t1443 * t1698 + 0.2e1 * t1447 * t1690 + 0.2e1 * t1451 * t1682, (-t1468 * t1647 - t1471 * t1646 - t1474 * t1645) * t1496, (-t1477 * t1647 - t1480 * t1646 - t1483 * t1645) * t1496, (t1347 ^ 2 / t1376 ^ 2 / 0.4e1 + t1346 ^ 2 / t1375 ^ 2 / 0.4e1 + t1345 ^ 2 / t1374 ^ 2 / 0.4e1) * t1497, 0, 0, t1401 * t1762 + t1402 * t1760 + t1403 * t1758, 0.4e1 * t1359 * t1762 + 0.4e1 * t1360 * t1760 + 0.4e1 * t1361 * t1758, -t1407 * t1669 - t1408 * t1666 - t1409 * t1664 + (t1409 * t1558 + t1408 * t1560 + t1407 * t1562 + (-t1362 * t1633 - t1363 * t1630 - t1364 * t1627) * t1495) * t1496, -t1404 * t1669 - t1405 * t1666 - t1406 * t1664 + (t1406 * t1558 + t1405 * t1560 + t1404 * t1562 + (-t1362 * t1642 - t1363 * t1639 - t1364 * t1636) * t1495) * t1496, (t1328 * t1661 + t1327 * t1662 + t1326 * t1663 + (t1326 * t1796 + t1327 * t1795 + t1328 * t1794) * t1495) * t1496, -t1239 * t1787 - t1241 * t1784 - t1243 * t1781 + (t1297 * t1661 + t1295 * t1662 + t1293 * t1663 + (t1314 * t1796 + t1315 * t1795 + t1316 * t1794) * t1495) * t1496, -t1240 * t1787 - t1242 * t1784 - t1244 * t1781 + (t1298 * t1661 + t1296 * t1662 + t1294 * t1663 + (t1317 * t1796 + t1318 * t1795 + t1319 * t1794) * t1495) * t1496, 1;];
tau_reg  = t1;
