% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRRRR7V1G3A0
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
% Datum: 2020-08-07 09:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V1G3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V1G3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 08:55:58
% EndTime: 2020-08-07 08:56:17
% DurationCPUTime: 20.40s
% Computational Cost: add. (18447->706), mult. (32262->1316), div. (3480->27), fcn. (26529->60), ass. (0->588)
t1449 = sin(qJ(1,3));
t1465 = -pkin(5) - pkin(4);
t1405 = t1449 * t1465;
t1458 = cos(qJ(1,3));
t1456 = cos(qJ(3,3));
t1796 = t1456 * pkin(2);
t1400 = pkin(1) + t1796;
t1457 = cos(qJ(2,3));
t1447 = sin(qJ(3,3));
t1448 = sin(qJ(2,3));
t1733 = t1447 * t1448;
t1525 = pkin(2) * t1733 - t1400 * t1457;
t1808 = t1458 * t1525 + t1405;
t1452 = sin(qJ(1,2));
t1406 = t1452 * t1465;
t1461 = cos(qJ(1,2));
t1459 = cos(qJ(3,2));
t1795 = t1459 * pkin(2);
t1402 = pkin(1) + t1795;
t1460 = cos(qJ(2,2));
t1450 = sin(qJ(3,2));
t1451 = sin(qJ(2,2));
t1729 = t1450 * t1451;
t1524 = pkin(2) * t1729 - t1402 * t1460;
t1807 = t1461 * t1524 + t1406;
t1455 = sin(qJ(1,1));
t1407 = t1455 * t1465;
t1464 = cos(qJ(1,1));
t1462 = cos(qJ(3,1));
t1794 = t1462 * pkin(2);
t1404 = pkin(1) + t1794;
t1463 = cos(qJ(2,1));
t1453 = sin(qJ(3,1));
t1454 = sin(qJ(2,1));
t1725 = t1453 * t1454;
t1523 = pkin(2) * t1725 - t1404 * t1463;
t1806 = t1464 * t1523 + t1407;
t1805 = -0.2e1 * pkin(1);
t1804 = 0.2e1 * pkin(1);
t1803 = 0.2e1 * pkin(2);
t1802 = 0.2e1 * t1465;
t1476 = 0.1e1 / pkin(1);
t1801 = 0.2e1 * t1476;
t1800 = pkin(4) / 0.2e1;
t1799 = pkin(4) * t1448;
t1798 = pkin(4) * t1451;
t1797 = pkin(4) * t1454;
t1793 = -qJ(3,1) + qJ(1,1);
t1792 = qJ(3,1) + qJ(1,1);
t1791 = -qJ(3,2) + qJ(1,2);
t1790 = qJ(3,2) + qJ(1,2);
t1789 = -qJ(3,3) + qJ(1,3);
t1788 = qJ(3,3) + qJ(1,3);
t1787 = qJ(1,1) + 0.2e1 * qJ(2,1);
t1786 = qJ(1,1) - 0.2e1 * qJ(2,1);
t1785 = qJ(1,2) + 0.2e1 * qJ(2,2);
t1784 = qJ(1,2) - 0.2e1 * qJ(2,2);
t1783 = 0.2e1 * qJ(2,3) + qJ(1,3);
t1782 = -0.2e1 * qJ(2,3) + qJ(1,3);
t1466 = 0.2e1 * qJ(3,3);
t1325 = (cos(qJ(2,3) - t1789) + cos(qJ(2,3) + t1788)) * t1802 + (-sin(0.2e1 * qJ(3,3) - t1782) + sin(t1466 + t1783) + 0.2e1 * t1449) * pkin(2) + (-sin(qJ(3,3) - t1782) + sin(qJ(3,3) + t1783) + sin(t1789) + sin(t1788)) * pkin(1);
t1441 = qJ(2,3) + qJ(3,3);
t1351 = (-sin(t1466 + qJ(2,3)) + t1448) * pkin(2) + (sin(qJ(2,3) - qJ(3,3)) - sin(t1441)) * pkin(1);
t1781 = t1325 / t1351;
t1469 = 0.2e1 * qJ(3,2);
t1326 = (cos(qJ(2,2) - t1791) + cos(qJ(2,2) + t1790)) * t1802 + (-sin(0.2e1 * qJ(3,2) - t1784) + sin(t1469 + t1785) + 0.2e1 * t1452) * pkin(2) + (-sin(qJ(3,2) - t1784) + sin(qJ(3,2) + t1785) + sin(t1791) + sin(t1790)) * pkin(1);
t1442 = qJ(2,2) + qJ(3,2);
t1352 = (-sin(t1469 + qJ(2,2)) + t1451) * pkin(2) + (sin(qJ(2,2) - qJ(3,2)) - sin(t1442)) * pkin(1);
t1780 = t1326 / t1352;
t1472 = 0.2e1 * qJ(3,1);
t1327 = (cos(qJ(2,1) - t1793) + cos(qJ(2,1) + t1792)) * t1802 + (-sin(0.2e1 * qJ(3,1) - t1786) + sin(t1472 + t1787) + 0.2e1 * t1455) * pkin(2) + (-sin(qJ(3,1) - t1786) + sin(qJ(3,1) + t1787) + sin(t1793) + sin(t1792)) * pkin(1);
t1443 = qJ(2,1) + qJ(3,1);
t1353 = (-sin(t1472 + qJ(2,1)) + t1454) * pkin(2) + (sin(qJ(2,1) - qJ(3,1)) - sin(t1443)) * pkin(1);
t1779 = t1327 / t1353;
t1732 = t1447 * t1457;
t1369 = pkin(2) * t1732 + t1448 * t1400;
t1444 = legFrame(3,2);
t1414 = sin(t1444);
t1417 = cos(t1444);
t1330 = -t1417 * t1369 - t1808 * t1414;
t1420 = 0.1e1 / t1447;
t1778 = t1330 * t1420;
t1331 = -t1414 * t1369 + t1808 * t1417;
t1777 = t1331 * t1420;
t1728 = t1450 * t1460;
t1370 = pkin(2) * t1728 + t1451 * t1402;
t1445 = legFrame(2,2);
t1415 = sin(t1445);
t1418 = cos(t1445);
t1332 = -t1418 * t1370 - t1807 * t1415;
t1424 = 0.1e1 / t1450;
t1776 = t1332 * t1424;
t1333 = -t1415 * t1370 + t1807 * t1418;
t1775 = t1333 * t1424;
t1724 = t1453 * t1463;
t1371 = pkin(2) * t1724 + t1454 * t1404;
t1446 = legFrame(1,2);
t1416 = sin(t1446);
t1419 = cos(t1446);
t1334 = -t1419 * t1371 - t1806 * t1416;
t1428 = 0.1e1 / t1453;
t1774 = t1334 * t1428;
t1335 = -t1416 * t1371 + t1806 * t1419;
t1773 = t1335 * t1428;
t1339 = -t1449 * t1525 + t1458 * t1465;
t1772 = t1339 * t1420;
t1340 = -t1452 * t1524 + t1461 * t1465;
t1771 = t1340 * t1424;
t1341 = -t1455 * t1523 + t1464 * t1465;
t1770 = t1341 * t1428;
t1363 = 0.1e1 / t1525;
t1769 = t1363 * t1420;
t1768 = 0.1e1 / t1525 ^ 2 / t1447 ^ 2;
t1365 = 0.1e1 / t1524;
t1767 = t1365 * t1424;
t1766 = 0.1e1 / t1524 ^ 2 / t1450 ^ 2;
t1367 = 0.1e1 / t1523;
t1765 = t1367 * t1428;
t1764 = 0.1e1 / t1523 ^ 2 / t1453 ^ 2;
t1393 = pkin(2) * cos(t1441) + t1457 * pkin(1);
t1387 = 0.1e1 / t1393;
t1763 = t1387 * t1449;
t1762 = t1387 * t1458;
t1388 = 0.1e1 / t1393 ^ 2;
t1423 = t1449 ^ 2;
t1761 = t1388 * t1423;
t1394 = pkin(2) * cos(t1442) + t1460 * pkin(1);
t1389 = 0.1e1 / t1394;
t1760 = t1389 * t1452;
t1759 = t1389 * t1461;
t1390 = 0.1e1 / t1394 ^ 2;
t1427 = t1452 ^ 2;
t1758 = t1390 * t1427;
t1395 = pkin(2) * cos(t1443) + t1463 * pkin(1);
t1391 = 0.1e1 / t1395;
t1757 = t1391 * t1455;
t1756 = t1391 * t1464;
t1392 = 0.1e1 / t1395 ^ 2;
t1431 = t1455 ^ 2;
t1755 = t1392 * t1431;
t1432 = t1456 ^ 2;
t1396 = pkin(1) * t1456 + t1432 * t1803 - pkin(2);
t1754 = t1396 * t1448;
t1753 = t1396 * t1458;
t1435 = t1459 ^ 2;
t1397 = pkin(1) * t1459 + t1435 * t1803 - pkin(2);
t1752 = t1397 * t1451;
t1751 = t1397 * t1461;
t1438 = t1462 ^ 2;
t1398 = t1462 * pkin(1) + t1438 * t1803 - pkin(2);
t1750 = t1398 * t1454;
t1749 = t1398 * t1464;
t1748 = (pkin(1) + 0.2e1 * t1796) * t1447;
t1747 = (pkin(1) + 0.2e1 * t1795) * t1450;
t1746 = (pkin(1) + 0.2e1 * t1794) * t1453;
t1745 = t1414 * t1449;
t1744 = t1415 * t1452;
t1743 = t1416 * t1455;
t1742 = t1417 * t1449;
t1741 = t1418 * t1452;
t1740 = t1419 * t1455;
t1739 = t1420 * t1458;
t1738 = t1424 * t1461;
t1737 = t1428 * t1464;
t1434 = t1458 ^ 2;
t1736 = t1434 * t1388;
t1437 = t1461 ^ 2;
t1735 = t1437 * t1390;
t1440 = t1464 ^ 2;
t1734 = t1440 * t1392;
t1731 = t1448 * t1457;
t1730 = t1448 * t1458;
t1727 = t1451 * t1460;
t1726 = t1451 * t1461;
t1723 = t1454 * t1463;
t1722 = t1454 * t1464;
t1721 = t1456 * t1447;
t1720 = t1456 * t1457;
t1719 = t1459 * t1450;
t1718 = t1459 * t1460;
t1717 = t1462 * t1453;
t1716 = t1462 * t1463;
t1475 = 0.1e1 / pkin(2);
t1715 = t1475 * t1476;
t1714 = -0.2e1 * t1732;
t1713 = -0.2e1 * t1728;
t1712 = -0.2e1 * t1724;
t1711 = -0.2e1 * t1720;
t1710 = -0.2e1 * t1718;
t1709 = -0.2e1 * t1716;
t1708 = pkin(2) * t1721;
t1707 = pkin(2) * t1719;
t1706 = pkin(2) * t1717;
t1705 = pkin(4) * t1387 * t1457;
t1704 = pkin(4) * t1389 * t1460;
t1703 = pkin(4) * t1391 * t1463;
t1355 = -t1406 * t1729 + (t1435 - 0.1e1) * t1461 * pkin(2);
t1436 = t1460 ^ 2;
t1648 = t1450 * t1726;
t1497 = pkin(1) * t1648 + (t1648 * t1803 + t1406) * t1459;
t1303 = (t1415 * t1747 + t1418 * t1751) * t1436 + (t1415 * t1752 - t1418 * t1497) * t1460 - t1355 * t1418 - t1415 * t1707;
t1689 = t1476 * t1767;
t1283 = t1303 * t1689;
t1651 = t1424 * t1715;
t1252 = t1333 * t1651 - t1283;
t1702 = t1252 * t1759;
t1354 = -t1405 * t1733 + (t1432 - 0.1e1) * t1458 * pkin(2);
t1433 = t1457 ^ 2;
t1649 = t1447 * t1730;
t1498 = pkin(1) * t1649 + (t1649 * t1803 + t1405) * t1456;
t1300 = (-t1414 * t1753 + t1417 * t1748) * t1433 + (t1414 * t1498 + t1417 * t1754) * t1457 + t1354 * t1414 - t1417 * t1708;
t1701 = t1300 * t1769;
t1301 = (t1414 * t1748 + t1417 * t1753) * t1433 + (t1414 * t1754 - t1417 * t1498) * t1457 - t1354 * t1417 - t1414 * t1708;
t1700 = t1301 * t1769;
t1302 = (-t1415 * t1751 + t1418 * t1747) * t1436 + (t1415 * t1497 + t1418 * t1752) * t1460 + t1355 * t1415 - t1418 * t1707;
t1699 = t1302 * t1767;
t1698 = t1303 * t1767;
t1356 = -t1407 * t1725 + (t1438 - 0.1e1) * t1464 * pkin(2);
t1439 = t1463 ^ 2;
t1647 = t1453 * t1722;
t1496 = pkin(1) * t1647 + (t1647 * t1803 + t1407) * t1462;
t1304 = (-t1416 * t1749 + t1419 * t1746) * t1439 + (t1416 * t1496 + t1419 * t1750) * t1463 + t1356 * t1416 - t1419 * t1706;
t1697 = t1304 * t1765;
t1305 = (t1416 * t1746 + t1419 * t1749) * t1439 + (t1416 * t1750 - t1419 * t1496) * t1463 - t1356 * t1419 - t1416 * t1706;
t1696 = t1305 * t1765;
t1643 = t1781 / 0.2e1;
t1315 = t1476 * t1643;
t1652 = t1420 * t1715;
t1306 = t1339 * t1652 + t1315;
t1695 = t1306 * t1769;
t1642 = t1780 / 0.2e1;
t1316 = t1476 * t1642;
t1307 = t1340 * t1651 + t1316;
t1694 = t1307 * t1767;
t1641 = t1779 / 0.2e1;
t1317 = t1476 * t1641;
t1650 = t1428 * t1715;
t1308 = t1341 * t1650 + t1317;
t1693 = t1308 * t1765;
t1692 = t1363 * t1739;
t1691 = t1476 * t1769;
t1690 = t1365 * t1738;
t1688 = t1367 * t1737;
t1687 = t1476 * t1765;
t1381 = t1720 - t1733;
t1686 = t1381 * t1762;
t1383 = t1716 - t1725;
t1685 = t1383 * t1756;
t1384 = t1456 * t1448 + t1732;
t1684 = t1384 * t1762;
t1386 = t1462 * t1454 + t1724;
t1683 = t1386 * t1756;
t1682 = t1387 * t1745;
t1681 = t1387 * t1742;
t1680 = t1420 * t1763;
t1679 = t1387 * t1739;
t1678 = t1448 * t1763;
t1677 = t1387 * t1730;
t1676 = t1388 * t1731;
t1675 = t1388 * t1449 * t1458;
t1674 = t1389 * t1744;
t1673 = t1389 * t1741;
t1672 = t1424 * t1760;
t1671 = t1389 * t1738;
t1670 = t1451 * t1760;
t1669 = t1389 * t1726;
t1668 = t1390 * t1727;
t1667 = t1390 * t1452 * t1461;
t1666 = t1391 * t1743;
t1665 = t1391 * t1740;
t1664 = t1428 * t1757;
t1663 = t1391 * t1737;
t1662 = t1454 * t1757;
t1661 = t1391 * t1722;
t1660 = t1392 * t1723;
t1659 = t1392 * t1455 * t1464;
t1408 = t1414 ^ 2;
t1658 = t1408 * t1761;
t1409 = t1415 ^ 2;
t1657 = t1409 * t1758;
t1410 = t1416 ^ 2;
t1656 = t1410 * t1755;
t1411 = t1417 ^ 2;
t1655 = t1411 * t1761;
t1412 = t1418 ^ 2;
t1654 = t1412 * t1758;
t1413 = t1419 ^ 2;
t1653 = t1413 * t1755;
t1280 = t1300 * t1691;
t1249 = t1330 * t1652 - t1280;
t1646 = t1249 * t1762;
t1282 = t1302 * t1689;
t1251 = t1332 * t1651 - t1282;
t1645 = t1251 * t1759;
t1644 = t1307 * t1759;
t1640 = t1715 / 0.2e1;
t1639 = t1387 * t1433 * t1804;
t1638 = t1389 * t1436 * t1804;
t1637 = t1391 * t1439 * t1804;
t1636 = t1449 * t1705;
t1635 = t1458 * t1705;
t1634 = t1452 * t1704;
t1633 = t1461 * t1704;
t1632 = t1455 * t1703;
t1631 = t1464 * t1703;
t1630 = pkin(4) * t1677;
t1629 = pkin(4) * t1669;
t1628 = pkin(4) * t1661;
t1627 = t1762 * t1781;
t1626 = t1759 * t1780;
t1625 = t1756 * t1779;
t1624 = t1381 * t1682;
t1623 = t1381 * t1681;
t1622 = t1381 * t1679;
t1382 = t1718 - t1729;
t1621 = t1382 * t1674;
t1620 = t1382 * t1673;
t1619 = t1382 * t1671;
t1618 = t1383 * t1666;
t1617 = t1383 * t1665;
t1616 = t1383 * t1663;
t1615 = t1384 * t1682;
t1614 = t1384 * t1681;
t1613 = t1384 * t1679;
t1385 = t1459 * t1451 + t1728;
t1612 = t1385 * t1674;
t1611 = t1385 * t1673;
t1610 = t1385 * t1671;
t1609 = t1386 * t1666;
t1608 = t1386 * t1665;
t1607 = t1386 * t1663;
t1606 = t1414 * t1680;
t1605 = t1414 * t1678;
t1604 = t1417 * t1680;
t1603 = t1417 * t1678;
t1602 = t1423 * t1676;
t1601 = t1415 * t1672;
t1600 = t1415 * t1670;
t1599 = t1418 * t1672;
t1598 = t1418 * t1670;
t1597 = t1427 * t1668;
t1596 = t1416 * t1664;
t1595 = t1416 * t1662;
t1594 = t1419 * t1664;
t1593 = t1419 * t1662;
t1592 = t1431 * t1660;
t1591 = t1414 * t1675;
t1590 = t1415 * t1667;
t1589 = t1416 * t1659;
t1588 = t1417 * t1414 * t1761;
t1587 = t1417 * t1675;
t1586 = t1418 * t1415 * t1758;
t1585 = t1418 * t1667;
t1584 = t1419 * t1416 * t1755;
t1583 = t1419 * t1659;
t1582 = t1420 * t1640;
t1581 = t1424 * t1640;
t1580 = t1428 * t1640;
t1579 = t1449 * t1639;
t1578 = t1452 * t1638;
t1577 = t1455 * t1637;
t1576 = t1447 * t1636;
t1575 = t1456 * t1636;
t1574 = t1450 * t1634;
t1573 = t1459 * t1634;
t1572 = t1453 * t1632;
t1571 = t1462 * t1632;
t1570 = pkin(4) * t1605;
t1569 = pkin(4) * t1603;
t1568 = pkin(4) * t1600;
t1567 = pkin(4) * t1598;
t1566 = pkin(4) * t1595;
t1565 = pkin(4) * t1593;
t1564 = t1363 * t1622;
t1563 = t1363 * t1613;
t1562 = t1365 * t1619;
t1561 = t1365 * t1610;
t1560 = t1367 * t1616;
t1559 = t1367 * t1607;
t1558 = t1381 * t1606;
t1557 = t1381 * t1604;
t1556 = t1382 * t1601;
t1555 = t1382 * t1599;
t1554 = t1383 * t1596;
t1553 = t1383 * t1594;
t1552 = t1384 * t1606;
t1551 = t1384 * t1604;
t1550 = t1385 * t1601;
t1549 = t1385 * t1599;
t1548 = t1386 * t1596;
t1547 = t1386 * t1594;
t1546 = t1675 * t1731;
t1545 = t1667 * t1727;
t1544 = t1659 * t1723;
t1543 = t1643 * t1769;
t1542 = -t1627 / 0.2e1;
t1541 = t1642 * t1767;
t1540 = -t1626 / 0.2e1;
t1539 = t1641 * t1765;
t1538 = -t1625 / 0.2e1;
t1537 = t1643 * t1745;
t1536 = t1642 * t1744;
t1535 = t1641 * t1743;
t1534 = -t1742 * t1781 / 0.2e1;
t1533 = -t1741 * t1780 / 0.2e1;
t1532 = -t1740 * t1779 / 0.2e1;
t1531 = t1414 * t1576;
t1530 = t1415 * t1574;
t1529 = t1416 * t1572;
t1528 = t1414 * t1575;
t1527 = t1415 * t1573;
t1526 = t1416 * t1571;
t1522 = t1300 * t1363 * t1606;
t1521 = t1301 * t1363 * t1604;
t1520 = t1302 * t1365 * t1601;
t1519 = t1303 * t1365 * t1599;
t1518 = t1304 * t1367 * t1596;
t1517 = t1305 * t1367 * t1594;
t1516 = t1363 * t1558;
t1515 = t1363 * t1557;
t1514 = t1363 * t1552;
t1513 = t1363 * t1551;
t1512 = t1365 * t1556;
t1511 = t1365 * t1555;
t1510 = t1365 * t1550;
t1509 = t1365 * t1549;
t1508 = t1367 * t1554;
t1507 = t1367 * t1553;
t1506 = t1367 * t1548;
t1505 = t1367 * t1547;
t1504 = t1387 * t1537;
t1503 = t1387 * t1534;
t1502 = t1389 * t1536;
t1501 = t1389 * t1533;
t1500 = t1391 * t1535;
t1499 = t1391 * t1532;
t1495 = t1306 * t1799 + t1458 * t1639;
t1494 = t1307 * t1798 + t1461 * t1638;
t1493 = t1308 * t1797 + t1464 * t1637;
t1492 = -t1249 * t1799 + t1414 * t1579;
t1281 = t1301 * t1691;
t1250 = t1331 * t1652 - t1281;
t1491 = t1250 * t1799 + t1417 * t1579;
t1490 = -t1251 * t1798 + t1415 * t1578;
t1489 = t1252 * t1798 + t1418 * t1578;
t1284 = t1304 * t1687;
t1253 = t1334 * t1650 - t1284;
t1488 = -t1253 * t1797 + t1416 * t1577;
t1285 = t1305 * t1687;
t1254 = t1335 * t1650 - t1285;
t1487 = t1254 * t1797 + t1419 * t1577;
t1486 = t1363 * (-t1300 * t1417 + t1301 * t1414) * t1680;
t1485 = t1365 * (-t1302 * t1418 + t1303 * t1415) * t1672;
t1484 = t1367 * (-t1304 * t1419 + t1305 * t1416) * t1664;
t1483 = t1387 * (t1300 * t1692 + t1537);
t1482 = t1387 * (t1301 * t1692 + t1534);
t1481 = t1389 * (t1302 * t1690 + t1536);
t1480 = t1389 * (t1303 * t1690 + t1533);
t1479 = t1391 * (t1304 * t1688 + t1535);
t1478 = t1391 * (t1305 * t1688 + t1532);
t1477 = 0.1e1 / pkin(1) ^ 2;
t1430 = t1454 ^ 2;
t1426 = t1451 ^ 2;
t1422 = t1448 ^ 2;
t1380 = t1386 ^ 2;
t1379 = t1385 ^ 2;
t1378 = t1384 ^ 2;
t1362 = t1462 * t1631;
t1361 = t1459 * t1633;
t1360 = t1456 * t1635;
t1359 = t1453 * t1631;
t1358 = t1450 * t1633;
t1357 = t1447 * t1635;
t1347 = t1419 * t1571;
t1346 = t1418 * t1573;
t1345 = t1417 * t1575;
t1344 = t1419 * t1572;
t1343 = t1418 * t1574;
t1342 = t1417 * t1576;
t1338 = (t1438 - 0.1e1 / 0.2e1) * t1723 + (t1439 - 0.1e1 / 0.2e1) * t1717;
t1337 = (t1435 - 0.1e1 / 0.2e1) * t1727 + (t1436 - 0.1e1 / 0.2e1) * t1719;
t1336 = (t1432 - 0.1e1 / 0.2e1) * t1731 + (t1433 - 0.1e1 / 0.2e1) * t1721;
t1329 = t1583 + t1585 + t1587;
t1328 = -t1589 - t1590 - t1591;
t1324 = -t1584 - t1586 - t1588;
t1323 = t1422 * t1587 + t1426 * t1585 + t1430 * t1583;
t1322 = -t1422 * t1591 - t1426 * t1590 - t1430 * t1589;
t1321 = -t1422 * t1588 - t1426 * t1586 - t1430 * t1584;
t1320 = 0.2e1 * t1417 * t1546 + 0.2e1 * t1418 * t1545 + 0.2e1 * t1419 * t1544;
t1319 = -0.2e1 * t1414 * t1546 - 0.2e1 * t1415 * t1545 - 0.2e1 * t1416 * t1544;
t1318 = -0.2e1 * t1584 * t1723 - 0.2e1 * t1586 * t1727 - 0.2e1 * t1588 * t1731;
t1314 = t1378 * t1587 + t1379 * t1585 + t1380 * t1583;
t1313 = -t1378 * t1591 - t1379 * t1590 - t1380 * t1589;
t1312 = -t1378 * t1588 - t1379 * t1586 - t1380 * t1584;
t1311 = t1628 + t1641;
t1310 = t1629 + t1642;
t1309 = t1630 + t1643;
t1299 = -t1453 * t1311 + t1362;
t1298 = -t1450 * t1310 + t1361;
t1297 = -t1447 * t1309 + t1360;
t1296 = t1462 * t1311 + t1359;
t1295 = t1459 * t1310 + t1358;
t1294 = t1456 * t1309 + t1357;
t1293 = -pkin(1) * t1661 + t1308 * t1800;
t1292 = t1628 + (t1341 * t1580 + t1317) * t1804;
t1291 = -pkin(1) * t1669 + t1307 * t1800;
t1290 = t1629 + (t1340 * t1581 + t1316) * t1804;
t1289 = -pkin(1) * t1677 + t1306 * t1800;
t1288 = t1630 + (t1339 * t1582 + t1315) * t1804;
t1287 = 0.4e1 * t1336 * t1587 + 0.4e1 * t1337 * t1585 + 0.4e1 * t1338 * t1583;
t1286 = -0.4e1 * t1336 * t1591 - 0.4e1 * t1337 * t1590 - 0.4e1 * t1338 * t1589;
t1279 = -0.4e1 * t1336 * t1588 - 0.4e1 * t1337 * t1586 - 0.4e1 * t1338 * t1584;
t1278 = -t1292 * t1453 + t1362;
t1277 = t1292 * t1462 + t1359;
t1276 = -t1290 * t1450 + t1361;
t1275 = t1290 * t1459 + t1358;
t1274 = -t1288 * t1447 + t1360;
t1273 = t1288 * t1456 + t1357;
t1272 = t1565 - t1696;
t1271 = -t1566 - t1697;
t1270 = t1567 - t1698;
t1269 = -t1568 - t1699;
t1268 = t1569 - t1700;
t1267 = -t1570 - t1701;
t1266 = -t1453 * t1272 + t1347;
t1265 = -t1453 * t1271 - t1526;
t1264 = -t1450 * t1270 + t1346;
t1263 = -t1450 * t1269 - t1527;
t1262 = -t1447 * t1268 + t1345;
t1261 = -t1447 * t1267 - t1528;
t1260 = t1462 * t1272 + t1344;
t1259 = t1462 * t1271 - t1529;
t1258 = t1459 * t1270 + t1343;
t1257 = t1459 * t1269 - t1530;
t1256 = t1456 * t1268 + t1342;
t1255 = t1456 * t1267 - t1531;
t1248 = -pkin(1) * t1593 + t1254 * t1800;
t1247 = pkin(1) * t1595 + t1253 * t1800;
t1246 = t1565 + (t1335 * t1580 - t1285) * t1804;
t1245 = t1566 + (t1334 * t1580 - t1284) * t1805;
t1244 = -pkin(1) * t1598 + t1252 * t1800;
t1243 = pkin(1) * t1600 + t1251 * t1800;
t1242 = t1567 + (t1333 * t1581 - t1283) * t1804;
t1241 = t1568 + (t1332 * t1581 - t1282) * t1805;
t1240 = -pkin(1) * t1603 + t1250 * t1800;
t1239 = pkin(1) * t1605 + t1249 * t1800;
t1238 = t1569 + (t1331 * t1582 - t1281) * t1804;
t1237 = t1570 + (t1330 * t1582 - t1280) * t1805;
t1236 = -t1246 * t1453 + t1347;
t1235 = t1246 * t1462 + t1344;
t1234 = -t1245 * t1462 - t1529;
t1233 = t1245 * t1453 - t1526;
t1232 = -t1242 * t1450 + t1346;
t1231 = t1242 * t1459 + t1343;
t1230 = -t1241 * t1459 - t1530;
t1229 = t1241 * t1450 - t1527;
t1228 = -t1238 * t1447 + t1345;
t1227 = t1238 * t1456 + t1342;
t1226 = -t1237 * t1456 - t1531;
t1225 = t1237 * t1447 - t1528;
t1224 = t1293 * t1709 + t1453 * t1493;
t1223 = t1293 * t1712 - t1462 * t1493;
t1222 = t1291 * t1710 + t1450 * t1494;
t1221 = t1291 * t1713 - t1459 * t1494;
t1220 = t1289 * t1711 + t1447 * t1495;
t1219 = t1289 * t1714 - t1456 * t1495;
t1218 = t1248 * t1709 + t1453 * t1487;
t1217 = t1247 * t1709 - t1453 * t1488;
t1216 = t1248 * t1712 - t1462 * t1487;
t1215 = t1247 * t1712 + t1462 * t1488;
t1214 = t1244 * t1710 + t1450 * t1489;
t1213 = t1243 * t1710 - t1450 * t1490;
t1212 = t1244 * t1713 - t1459 * t1489;
t1211 = t1243 * t1713 + t1459 * t1490;
t1210 = t1240 * t1711 + t1447 * t1491;
t1209 = t1239 * t1711 - t1447 * t1492;
t1208 = t1240 * t1714 - t1456 * t1491;
t1207 = t1239 * t1714 + t1456 * t1492;
t1206 = (-t1301 * t1543 - t1303 * t1541 - t1305 * t1539) * t1477;
t1205 = (-t1300 * t1543 - t1302 * t1541 - t1304 * t1539) * t1477;
t1204 = (t1300 * t1301 * t1768 + t1302 * t1303 * t1766 + t1304 * t1305 * t1764) * t1477;
t1203 = (t1457 * t1482 + t1460 * t1480 + t1463 * t1478) * t1476;
t1202 = (t1448 * t1482 + t1451 * t1480 + t1454 * t1478) * t1476;
t1201 = (t1457 * t1483 + t1460 * t1481 + t1463 * t1479) * t1476;
t1200 = (t1448 * t1483 + t1451 * t1481 + t1454 * t1479) * t1476;
t1199 = (-t1457 * t1486 - t1460 * t1485 - t1463 * t1484) * t1476;
t1198 = (-t1448 * t1486 - t1451 * t1485 - t1454 * t1484) * t1476;
t1 = [t1653 + t1654 + t1655, 0, 0, t1422 * t1655 + t1426 * t1654 + t1430 * t1653, 0.2e1 * t1411 * t1602 + 0.2e1 * t1412 * t1597 + 0.2e1 * t1413 * t1592, (t1448 * t1521 + t1451 * t1519 + t1454 * t1517) * t1801, (t1457 * t1521 + t1460 * t1519 + t1463 * t1517) * t1801, (t1301 ^ 2 * t1768 + t1303 ^ 2 * t1766 + t1305 ^ 2 * t1764) * t1477, 0, 0, t1378 * t1655 + t1379 * t1654 + t1380 * t1653, 0.4e1 * t1336 * t1655 + 0.4e1 * t1337 * t1654 + 0.4e1 * t1338 * t1653, -t1250 * t1614 - t1252 * t1611 - t1254 * t1608 + (t1301 * t1513 + t1303 * t1509 + t1305 * t1505 + (-t1331 * t1551 - t1333 * t1549 - t1335 * t1547) * t1475) * t1476, -t1250 * t1623 - t1252 * t1620 - t1254 * t1617 + (t1301 * t1515 + t1303 * t1511 + t1305 * t1507 + (-t1331 * t1557 - t1333 * t1555 - t1335 * t1553) * t1475) * t1476, (-t1250 * t1700 - t1252 * t1698 - t1254 * t1696 + (t1250 * t1777 + t1252 * t1775 + t1254 * t1773) * t1475) * t1476, -t1208 * t1681 - t1212 * t1673 - t1216 * t1665 + (-t1227 * t1700 - t1231 * t1698 - t1235 * t1696 + (t1256 * t1777 + t1258 * t1775 + t1260 * t1773) * t1475) * t1476, -t1210 * t1681 - t1214 * t1673 - t1218 * t1665 + (-t1228 * t1700 - t1232 * t1698 - t1236 * t1696 + (t1262 * t1777 + t1264 * t1775 + t1266 * t1773) * t1475) * t1476, 1; t1324, 0, 0, t1321, t1318, t1198, t1199, t1204, 0, 0, t1312, t1279, -t1249 * t1614 - t1251 * t1611 - t1253 * t1608 + (-t1301 * t1514 - t1303 * t1510 - t1305 * t1506 + (t1331 * t1552 + t1333 * t1550 + t1335 * t1548) * t1475) * t1476, -t1249 * t1623 - t1251 * t1620 - t1253 * t1617 + (-t1301 * t1516 - t1303 * t1512 - t1305 * t1508 + (t1331 * t1558 + t1333 * t1556 + t1335 * t1554) * t1475) * t1476, (-t1249 * t1700 - t1251 * t1698 - t1253 * t1696 + (t1249 * t1777 + t1251 * t1775 + t1253 * t1773) * t1475) * t1476, -t1207 * t1681 - t1211 * t1673 - t1215 * t1665 + (-t1226 * t1700 - t1230 * t1698 - t1234 * t1696 + (t1255 * t1777 + t1257 * t1775 + t1259 * t1773) * t1475) * t1476, -t1209 * t1681 - t1213 * t1673 - t1217 * t1665 + (-t1225 * t1700 - t1229 * t1698 - t1233 * t1696 + (t1261 * t1777 + t1263 * t1775 + t1265 * t1773) * t1475) * t1476, 0; t1329, 0, 0, t1323, t1320, t1202, t1203, t1206, 0, 0, t1314, t1287, -t1306 * t1614 - t1307 * t1611 - t1308 * t1608 + (t1301 * t1563 + t1303 * t1561 + t1305 * t1559 + (-t1331 * t1613 - t1333 * t1610 - t1335 * t1607) * t1475) * t1476, -t1306 * t1623 - t1307 * t1620 - t1308 * t1617 + (t1301 * t1564 + t1303 * t1562 + t1305 * t1560 + (-t1331 * t1622 - t1333 * t1619 - t1335 * t1616) * t1475) * t1476, (-t1301 * t1695 - t1303 * t1694 - t1305 * t1693 + (t1306 * t1777 + t1307 * t1775 + t1308 * t1773) * t1475) * t1476, -t1219 * t1681 - t1221 * t1673 - t1223 * t1665 + (-t1273 * t1700 - t1275 * t1698 - t1277 * t1696 + (t1294 * t1777 + t1295 * t1775 + t1296 * t1773) * t1475) * t1476, -t1220 * t1681 - t1222 * t1673 - t1224 * t1665 + (-t1274 * t1700 - t1276 * t1698 - t1278 * t1696 + (t1297 * t1777 + t1298 * t1775 + t1299 * t1773) * t1475) * t1476, 0; t1324, 0, 0, t1321, t1318, t1198, t1199, t1204, 0, 0, t1312, t1279, t1250 * t1615 + t1252 * t1612 + t1254 * t1609 + (t1300 * t1513 + t1302 * t1509 + t1304 * t1505 + (-t1330 * t1551 - t1332 * t1549 - t1334 * t1547) * t1475) * t1476, t1250 * t1624 + t1252 * t1621 + t1254 * t1618 + (t1300 * t1515 + t1302 * t1511 + t1304 * t1507 + (-t1330 * t1557 - t1332 * t1555 - t1334 * t1553) * t1475) * t1476, (-t1250 * t1701 - t1252 * t1699 - t1254 * t1697 + (t1250 * t1778 + t1252 * t1776 + t1254 * t1774) * t1475) * t1476, t1208 * t1682 + t1212 * t1674 + t1216 * t1666 + (-t1227 * t1701 - t1231 * t1699 - t1235 * t1697 + (t1256 * t1778 + t1258 * t1776 + t1260 * t1774) * t1475) * t1476, t1210 * t1682 + t1214 * t1674 + t1218 * t1666 + (-t1228 * t1701 - t1232 * t1699 - t1236 * t1697 + (t1262 * t1778 + t1264 * t1776 + t1266 * t1774) * t1475) * t1476, 0; t1656 + t1657 + t1658, 0, 0, t1422 * t1658 + t1426 * t1657 + t1430 * t1656, 0.2e1 * t1408 * t1602 + 0.2e1 * t1409 * t1597 + 0.2e1 * t1410 * t1592, (-t1448 * t1522 - t1451 * t1520 - t1454 * t1518) * t1801, (-t1457 * t1522 - t1460 * t1520 - t1463 * t1518) * t1801, (t1300 ^ 2 * t1768 + t1302 ^ 2 * t1766 + t1304 ^ 2 * t1764) * t1477, 0, 0, t1378 * t1658 + t1379 * t1657 + t1380 * t1656, 0.4e1 * t1336 * t1658 + 0.4e1 * t1337 * t1657 + 0.4e1 * t1338 * t1656, t1249 * t1615 + t1251 * t1612 + t1253 * t1609 + (-t1300 * t1514 - t1302 * t1510 - t1304 * t1506 + (t1330 * t1552 + t1332 * t1550 + t1334 * t1548) * t1475) * t1476, t1249 * t1624 + t1251 * t1621 + t1253 * t1618 + (-t1300 * t1516 - t1302 * t1512 - t1304 * t1508 + (t1330 * t1558 + t1332 * t1556 + t1334 * t1554) * t1475) * t1476, (-t1249 * t1701 - t1251 * t1699 - t1253 * t1697 + (t1249 * t1778 + t1251 * t1776 + t1253 * t1774) * t1475) * t1476, t1207 * t1682 + t1211 * t1674 + t1215 * t1666 + (-t1226 * t1701 - t1230 * t1699 - t1234 * t1697 + (t1255 * t1778 + t1257 * t1776 + t1259 * t1774) * t1475) * t1476, t1209 * t1682 + t1213 * t1674 + t1217 * t1666 + (-t1225 * t1701 - t1229 * t1699 - t1233 * t1697 + (t1261 * t1778 + t1263 * t1776 + t1265 * t1774) * t1475) * t1476, 1; t1328, 0, 0, t1322, t1319, t1200, t1201, t1205, 0, 0, t1313, t1286, t1306 * t1615 + t1307 * t1612 + t1308 * t1609 + (t1300 * t1563 + t1302 * t1561 + t1304 * t1559 + (-t1330 * t1613 - t1332 * t1610 - t1334 * t1607) * t1475) * t1476, t1306 * t1624 + t1307 * t1621 + t1308 * t1618 + (t1300 * t1564 + t1302 * t1562 + t1304 * t1560 + (-t1330 * t1622 - t1332 * t1619 - t1334 * t1616) * t1475) * t1476, (-t1300 * t1695 - t1302 * t1694 - t1304 * t1693 + (t1306 * t1778 + t1307 * t1776 + t1308 * t1774) * t1475) * t1476, t1219 * t1682 + t1221 * t1674 + t1223 * t1666 + (-t1273 * t1701 - t1275 * t1699 - t1277 * t1697 + (t1294 * t1778 + t1295 * t1776 + t1296 * t1774) * t1475) * t1476, t1220 * t1682 + t1222 * t1674 + t1224 * t1666 + (-t1274 * t1701 - t1276 * t1699 - t1278 * t1697 + (t1297 * t1778 + t1298 * t1776 + t1299 * t1774) * t1475) * t1476, 0; t1329, 0, 0, t1323, t1320, t1202, t1203, t1206, 0, 0, t1314, t1287, -t1250 * t1684 - t1385 * t1702 - t1254 * t1683 + (t1386 * t1499 + t1385 * t1501 + t1384 * t1503 + (-t1339 * t1551 - t1340 * t1549 - t1341 * t1547) * t1475) * t1476, -t1250 * t1686 - t1382 * t1702 - t1254 * t1685 + (t1383 * t1499 + t1382 * t1501 + t1381 * t1503 + (-t1339 * t1557 - t1340 * t1555 - t1341 * t1553) * t1475) * t1476, (t1254 * t1641 + t1252 * t1642 + t1250 * t1643 + (t1250 * t1772 + t1252 * t1771 + t1254 * t1770) * t1475) * t1476, -t1208 * t1762 - t1212 * t1759 - t1216 * t1756 + (t1235 * t1641 + t1231 * t1642 + t1227 * t1643 + (t1256 * t1772 + t1258 * t1771 + t1260 * t1770) * t1475) * t1476, -t1210 * t1762 - t1214 * t1759 - t1218 * t1756 + (t1236 * t1641 + t1232 * t1642 + t1228 * t1643 + (t1262 * t1772 + t1264 * t1771 + t1266 * t1770) * t1475) * t1476, 0; t1328, 0, 0, t1322, t1319, t1200, t1201, t1205, 0, 0, t1313, t1286, -t1384 * t1646 - t1385 * t1645 - t1253 * t1683 + (t1386 * t1500 + t1385 * t1502 + t1384 * t1504 + (t1339 * t1552 + t1340 * t1550 + t1341 * t1548) * t1475) * t1476, -t1381 * t1646 - t1382 * t1645 - t1253 * t1685 + (t1383 * t1500 + t1382 * t1502 + t1381 * t1504 + (t1339 * t1558 + t1340 * t1556 + t1341 * t1554) * t1475) * t1476, (t1253 * t1641 + t1251 * t1642 + t1249 * t1643 + (t1249 * t1772 + t1251 * t1771 + t1253 * t1770) * t1475) * t1476, -t1207 * t1762 - t1211 * t1759 - t1215 * t1756 + (t1234 * t1641 + t1230 * t1642 + t1226 * t1643 + (t1255 * t1772 + t1257 * t1771 + t1259 * t1770) * t1475) * t1476, -t1209 * t1762 - t1213 * t1759 - t1217 * t1756 + (t1233 * t1641 + t1229 * t1642 + t1225 * t1643 + (t1261 * t1772 + t1263 * t1771 + t1265 * t1770) * t1475) * t1476, 0; t1734 + t1735 + t1736, 0, 0, t1422 * t1736 + t1426 * t1735 + t1430 * t1734, 0.2e1 * t1434 * t1676 + 0.2e1 * t1437 * t1668 + 0.2e1 * t1440 * t1660, (-t1448 * t1627 - t1451 * t1626 - t1454 * t1625) * t1476, (-t1457 * t1627 - t1460 * t1626 - t1463 * t1625) * t1476, (t1327 ^ 2 / t1353 ^ 2 / 0.4e1 + t1326 ^ 2 / t1352 ^ 2 / 0.4e1 + t1325 ^ 2 / t1351 ^ 2 / 0.4e1) * t1477, 0, 0, t1378 * t1736 + t1379 * t1735 + t1380 * t1734, 0.4e1 * t1336 * t1736 + 0.4e1 * t1337 * t1735 + 0.4e1 * t1338 * t1734, -t1306 * t1684 - t1385 * t1644 - t1308 * t1683 + (t1386 * t1538 + t1385 * t1540 + t1384 * t1542 + (-t1339 * t1613 - t1340 * t1610 - t1341 * t1607) * t1475) * t1476, -t1306 * t1686 - t1382 * t1644 - t1308 * t1685 + (t1383 * t1538 + t1382 * t1540 + t1381 * t1542 + (-t1339 * t1622 - t1340 * t1619 - t1341 * t1616) * t1475) * t1476, (t1308 * t1641 + t1307 * t1642 + t1306 * t1643 + (t1306 * t1772 + t1307 * t1771 + t1308 * t1770) * t1475) * t1476, -t1219 * t1762 - t1221 * t1759 - t1223 * t1756 + (t1277 * t1641 + t1275 * t1642 + t1273 * t1643 + (t1294 * t1772 + t1295 * t1771 + t1296 * t1770) * t1475) * t1476, -t1220 * t1762 - t1222 * t1759 - t1224 * t1756 + (t1278 * t1641 + t1276 * t1642 + t1274 * t1643 + (t1297 * t1772 + t1298 * t1771 + t1299 * t1770) * t1475) * t1476, 1;];
tau_reg  = t1;