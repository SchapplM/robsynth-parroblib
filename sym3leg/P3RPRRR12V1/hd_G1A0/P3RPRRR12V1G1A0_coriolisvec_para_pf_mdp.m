% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR12V1G1A0
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [14x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR12V1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(14,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'P3RPRRR12V1G1A0_coriolisvec_para_pf_mdp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:43
% EndTime: 2020-08-06 18:21:47
% DurationCPUTime: 3.59s
% Computational Cost: add. (10359->286), mult. (16168->529), div. (1761->18), fcn. (11982->18), ass. (0->231)
t1564 = sin(qJ(3,3));
t1523 = t1564 * pkin(3) + qJ(2,3);
t1513 = 0.1e1 / t1523 ^ 2;
t1570 = cos(qJ(3,3));
t1576 = xDP(3);
t1672 = t1570 * t1576;
t1561 = legFrame(3,3);
t1530 = sin(t1561);
t1533 = cos(t1561);
t1565 = sin(qJ(1,3));
t1571 = cos(qJ(1,3));
t1577 = xDP(2);
t1578 = xDP(1);
t1621 = t1565 * t1577 + t1571 * t1578;
t1622 = -t1565 * t1578 + t1571 * t1577;
t1491 = t1622 * t1530 + t1621 * t1533;
t1545 = 0.1e1 / t1564;
t1703 = t1491 * t1545;
t1579 = pkin(1) + pkin(5);
t1544 = pkin(6) + t1579;
t1704 = t1491 * t1544;
t1512 = 0.1e1 / t1523;
t1705 = t1491 * t1512;
t1521 = t1544 * t1578;
t1506 = -qJ(2,3) * t1577 + t1521;
t1522 = t1577 * t1544;
t1509 = qJ(2,3) * t1578 + t1522;
t1470 = ((t1506 * t1571 + t1565 * t1509) * t1533 + (-t1565 * t1506 + t1509 * t1571) * t1530) * t1564 + qJ(2,3) * t1672 + (t1564 * t1672 + (-t1621 * t1530 + t1622 * t1533) * (t1570 - 0.1e1) * (t1570 + 0.1e1)) * pkin(3);
t1718 = t1470 * t1545;
t1514 = t1512 * t1513;
t1719 = t1470 * t1514;
t1455 = (0.2e1 * t1513 * t1672 - t1719) * t1703 - (-t1704 + t1718) * t1513 * t1705;
t1595 = t1564 ^ 2;
t1546 = 0.1e1 / t1595;
t1560 = t1576 ^ 2;
t1585 = 0.1e1 / pkin(3) ^ 2;
t1676 = t1560 * t1585;
t1636 = t1579 * t1676;
t1629 = t1546 * t1636;
t1693 = t1513 * t1545;
t1725 = 0.2e1 * qJ(2,3);
t1737 = 0.2e1 * t1491 * t1470 * t1693 + t1455 * t1725 + t1629;
t1566 = sin(qJ(3,2));
t1524 = t1566 * pkin(3) + qJ(2,2);
t1516 = 0.1e1 / t1524 ^ 2;
t1572 = cos(qJ(3,2));
t1671 = t1572 * t1576;
t1562 = legFrame(2,3);
t1531 = sin(t1562);
t1534 = cos(t1562);
t1567 = sin(qJ(1,2));
t1573 = cos(qJ(1,2));
t1619 = t1567 * t1577 + t1573 * t1578;
t1620 = -t1567 * t1578 + t1573 * t1577;
t1492 = t1620 * t1531 + t1619 * t1534;
t1549 = 0.1e1 / t1566;
t1700 = t1492 * t1549;
t1701 = t1492 * t1544;
t1515 = 0.1e1 / t1524;
t1702 = t1492 * t1515;
t1507 = -qJ(2,2) * t1577 + t1521;
t1510 = qJ(2,2) * t1578 + t1522;
t1471 = ((t1507 * t1573 + t1567 * t1510) * t1534 + (-t1567 * t1507 + t1510 * t1573) * t1531) * t1566 + qJ(2,2) * t1671 + (t1566 * t1671 + (-t1619 * t1531 + t1620 * t1534) * (t1572 - 0.1e1) * (t1572 + 0.1e1)) * pkin(3);
t1716 = t1471 * t1549;
t1517 = t1515 * t1516;
t1717 = t1471 * t1517;
t1456 = (0.2e1 * t1516 * t1671 - t1717) * t1700 - (-t1701 + t1716) * t1516 * t1702;
t1598 = t1566 ^ 2;
t1550 = 0.1e1 / t1598;
t1628 = t1550 * t1636;
t1689 = t1516 * t1549;
t1726 = 0.2e1 * qJ(2,2);
t1736 = 0.2e1 * t1492 * t1471 * t1689 + t1456 * t1726 + t1628;
t1568 = sin(qJ(3,1));
t1525 = -t1568 * pkin(3) - qJ(2,1);
t1519 = 0.1e1 / t1525 ^ 2;
t1574 = cos(qJ(3,1));
t1670 = t1574 * t1576;
t1563 = legFrame(1,3);
t1532 = sin(t1563);
t1535 = cos(t1563);
t1569 = sin(qJ(1,1));
t1575 = cos(qJ(1,1));
t1617 = t1569 * t1577 + t1575 * t1578;
t1618 = -t1569 * t1578 + t1575 * t1577;
t1493 = t1618 * t1532 + t1617 * t1535;
t1553 = 0.1e1 / t1568;
t1697 = t1493 * t1553;
t1698 = t1493 * t1544;
t1518 = 0.1e1 / t1525;
t1699 = t1493 * t1518;
t1508 = -qJ(2,1) * t1577 + t1521;
t1511 = qJ(2,1) * t1578 + t1522;
t1472 = ((t1508 * t1575 + t1569 * t1511) * t1535 + (-t1569 * t1508 + t1511 * t1575) * t1532) * t1568 + qJ(2,1) * t1670 + (t1568 * t1670 + (-t1617 * t1532 + t1618 * t1535) * (t1574 - 0.1e1) * (t1574 + 0.1e1)) * pkin(3);
t1714 = t1472 * t1553;
t1520 = t1518 * t1519;
t1715 = t1472 * t1520;
t1457 = (0.2e1 * t1519 * t1670 + t1715) * t1697 + (-t1698 + t1714) * t1519 * t1699;
t1601 = t1568 ^ 2;
t1554 = 0.1e1 / t1601;
t1627 = t1554 * t1636;
t1686 = t1519 * t1553;
t1727 = 0.2e1 * qJ(2,1);
t1735 = 0.2e1 * t1493 * t1472 * t1686 + t1457 * t1727 + t1627;
t1557 = t1570 ^ 2;
t1734 = t1545 * t1557;
t1558 = t1572 ^ 2;
t1733 = t1549 * t1558;
t1559 = t1574 ^ 2;
t1732 = t1553 * t1559;
t1678 = t1554 * t1732;
t1637 = t1518 * t1678;
t1731 = t1518 * t1553 + t1637;
t1681 = t1550 * t1733;
t1639 = t1515 * t1681;
t1730 = -t1515 * t1549 - t1639;
t1684 = t1546 * t1734;
t1641 = t1512 * t1684;
t1729 = -t1512 * t1545 - t1641;
t1580 = qJ(2,3) ^ 2;
t1583 = pkin(3) ^ 2;
t1584 = 0.1e1 / pkin(3);
t1669 = t1576 * t1584;
t1623 = pkin(3) * t1544 * t1669;
t1658 = qJ(2,3) * t1669;
t1668 = -t1544 ^ 2 - t1583;
t1675 = t1564 * t1570;
t1694 = t1512 * t1570;
t1613 = -((-t1545 * t1623 * t1675 + (t1544 * t1564 * t1718 + ((t1557 * t1583 - t1580 + t1668) * t1564 + (t1557 - 0.1e1) * t1725 * pkin(3)) * t1491) * t1512) * t1513 + t1544 * t1719) * t1703 + ((t1545 * t1576 + t1694 * t1704) * t1564 + t1545 * t1658) * t1512 * t1546 * t1576;
t1581 = qJ(2,2) ^ 2;
t1660 = qJ(2,2) * t1669;
t1674 = t1566 * t1572;
t1690 = t1515 * t1572;
t1612 = -((-t1549 * t1623 * t1674 + (t1544 * t1566 * t1716 + ((t1558 * t1583 - t1581 + t1668) * t1566 + (t1558 - 0.1e1) * t1726 * pkin(3)) * t1492) * t1515) * t1516 + t1544 * t1717) * t1700 + ((t1549 * t1576 + t1690 * t1701) * t1566 + t1549 * t1660) * t1515 * t1550 * t1576;
t1582 = qJ(2,1) ^ 2;
t1662 = qJ(2,1) * t1669;
t1673 = t1568 * t1574;
t1687 = t1518 * t1574;
t1611 = (-(-t1553 * t1623 * t1673 - (t1544 * t1568 * t1714 + ((t1559 * t1583 - t1582 + t1668) * t1568 + (t1559 - 0.1e1) * t1727 * pkin(3)) * t1493) * t1518) * t1519 + t1544 * t1715) * t1697 - ((t1553 * t1576 - t1687 * t1698) * t1568 + t1553 * t1662) * t1518 * t1554 * t1576;
t1728 = -2 * MDP(8);
t1724 = pkin(1) * t1455;
t1723 = pkin(1) * t1456;
t1722 = pkin(1) * t1457;
t1721 = t1455 * t1512;
t1720 = t1456 * t1515;
t1496 = -t1569 * t1525 + t1544 * t1575;
t1499 = t1525 * t1575 + t1544 * t1569;
t1484 = t1496 * t1535 - t1532 * t1499;
t1713 = t1484 * t1518;
t1487 = t1532 * t1496 + t1499 * t1535;
t1712 = t1487 * t1518;
t1488 = t1491 ^ 2;
t1711 = t1488 * t1513;
t1710 = t1488 * t1514;
t1489 = t1492 ^ 2;
t1709 = t1489 * t1516;
t1708 = t1489 * t1517;
t1490 = t1493 ^ 2;
t1707 = t1490 * t1519;
t1706 = t1490 * t1520;
t1443 = t1613 - t1724;
t1696 = t1512 * t1443;
t1445 = t1612 - t1723;
t1692 = t1515 * t1445;
t1685 = t1545 * t1570;
t1683 = 0.1e1 / t1595 ^ 2 * t1570;
t1682 = t1549 * t1572;
t1680 = 0.1e1 / t1598 ^ 2 * t1572;
t1679 = t1553 * t1574;
t1677 = 0.1e1 / t1601 ^ 2 * t1574;
t1664 = 0.2e1 * t1669;
t1663 = qJ(2,1) * t1707;
t1661 = qJ(2,2) * t1709;
t1659 = qJ(2,3) * t1711;
t1479 = -t1546 * t1676 - t1711;
t1657 = t1479 * t1512 * t1564;
t1656 = t1479 * t1694;
t1480 = -t1550 * t1676 - t1709;
t1655 = t1480 * t1515 * t1566;
t1654 = t1480 * t1690;
t1481 = -t1554 * t1676 - t1707;
t1653 = t1481 * t1518 * t1568;
t1652 = t1481 * t1687;
t1494 = t1565 * t1523 + t1544 * t1571;
t1497 = -t1523 * t1571 + t1544 * t1565;
t1482 = t1494 * t1533 - t1530 * t1497;
t1651 = t1482 * t1710;
t1495 = t1567 * t1524 + t1544 * t1573;
t1498 = -t1524 * t1573 + t1544 * t1567;
t1483 = t1495 * t1534 - t1531 * t1498;
t1650 = t1483 * t1708;
t1649 = t1484 * t1706;
t1485 = t1530 * t1494 + t1497 * t1533;
t1648 = t1485 * t1710;
t1486 = t1531 * t1495 + t1498 * t1534;
t1647 = t1486 * t1708;
t1646 = t1487 * t1706;
t1645 = t1570 * t1711;
t1644 = t1572 * t1709;
t1643 = t1574 * t1707;
t1642 = t1546 * t1694;
t1640 = t1550 * t1690;
t1638 = t1554 * t1687;
t1635 = t1513 * t1664;
t1634 = t1516 * t1664;
t1633 = t1519 * t1664;
t1632 = t1545 * t1645;
t1631 = t1549 * t1644;
t1630 = t1553 * t1643;
t1626 = t1658 * t1705;
t1625 = t1660 * t1702;
t1624 = t1662 * t1699;
t1610 = t1455 * t1685 + t1456 * t1682 + t1457 * t1679;
t1527 = 0.2e1 * t1557 - 0.1e1;
t1587 = pkin(1) ^ 2;
t1609 = (t1570 * MDP(7) * t1635 + (MDP(8) * t1527 * t1635 + 0.2e1 * (MDP(6) * qJ(2,3) + MDP(5)) * t1719) * t1545) * t1491 + ((t1613 - 0.2e1 * t1724) * MDP(4) + ((t1580 + t1587) * t1455 - pkin(1) * t1613) * MDP(6) + (t1737 * t1564 - 0.2e1 * t1626 * t1685 + t1636 * t1684) * MDP(12) + (0.2e1 * t1626 + (-t1629 + t1737) * t1570) * MDP(13) + (MDP(5) * t1725 + t1557 * MDP(7) + t1675 * t1728 + MDP(1)) * t1455) * t1512;
t1528 = 0.2e1 * t1558 - 0.1e1;
t1608 = (t1572 * MDP(7) * t1634 + (MDP(8) * t1528 * t1634 + 0.2e1 * (MDP(6) * qJ(2,2) + MDP(5)) * t1717) * t1549) * t1492 + ((t1612 - 0.2e1 * t1723) * MDP(4) + ((t1581 + t1587) * t1456 - pkin(1) * t1612) * MDP(6) + (t1736 * t1566 - 0.2e1 * t1625 * t1682 + t1636 * t1681) * MDP(12) + (0.2e1 * t1625 + (-t1628 + t1736) * t1572) * MDP(13) + (MDP(5) * t1726 + t1558 * MDP(7) + t1674 * t1728 + MDP(1)) * t1456) * t1515;
t1529 = 0.2e1 * t1559 - 0.1e1;
t1607 = (t1574 * MDP(7) * t1633 + (MDP(8) * t1529 * t1633 - 0.2e1 * (MDP(6) * qJ(2,1) + MDP(5)) * t1715) * t1553) * t1493 - ((t1611 - 0.2e1 * t1722) * MDP(4) + ((t1582 + t1587) * t1457 - pkin(1) * t1611) * MDP(6) + (t1735 * t1568 + 0.2e1 * t1624 * t1679 + t1636 * t1678) * MDP(12) + (-0.2e1 * t1624 + (-t1627 + t1735) * t1574) * MDP(13) + (MDP(5) * t1727 + t1559 * MDP(7) + t1673 * t1728 + MDP(1)) * t1457) * t1518;
t1505 = t1532 * t1575 + t1535 * t1569;
t1504 = t1531 * t1573 + t1534 * t1567;
t1503 = t1530 * t1571 + t1533 * t1565;
t1502 = -t1532 * t1569 + t1535 * t1575;
t1501 = -t1531 * t1567 + t1534 * t1573;
t1500 = -t1530 * t1565 + t1533 * t1571;
t1447 = t1611 - t1722;
t1442 = t1579 * t1457 - t1611;
t1441 = t1579 * t1456 - t1612;
t1440 = t1579 * t1455 - t1613;
t1 = [(-t1457 * t1713 + t1482 * t1721 + t1483 * t1720) * MDP(4) + (t1649 - t1650 - t1651) * MDP(5) + (qJ(2,1) * t1649 - qJ(2,2) * t1650 - qJ(2,3) * t1651 - t1447 * t1713 + t1482 * t1696 + t1483 * t1692) * MDP(6) + (t1482 * t1657 + t1483 * t1655 - t1484 * t1653) * MDP(12) + (t1482 * t1656 + t1483 * t1654 - t1484 * t1652) * MDP(13) + t1607 * t1502 + t1608 * t1501 + t1609 * t1500 + ((t1729 * t1500 + t1730 * t1501 + t1731 * t1502) * MDP(9) + (-t1482 * t1641 - t1483 * t1639 + t1484 * t1637) * MDP(12) + (t1482 * t1642 + t1483 * t1640 - t1484 * t1638) * MDP(13)) * t1676; (-t1457 * t1712 + t1485 * t1721 + t1486 * t1720) * MDP(4) + (t1646 - t1647 - t1648) * MDP(5) + (qJ(2,1) * t1646 - qJ(2,2) * t1647 - qJ(2,3) * t1648 - t1447 * t1712 + t1485 * t1696 + t1486 * t1692) * MDP(6) + (t1485 * t1657 + t1486 * t1655 - t1487 * t1653) * MDP(12) + (t1485 * t1656 + t1486 * t1654 - t1487 * t1652) * MDP(13) + t1607 * t1505 + t1608 * t1504 + t1609 * t1503 + ((t1729 * t1503 + t1730 * t1504 + t1731 * t1505) * MDP(9) + (-t1485 * t1641 - t1486 * t1639 + t1487 * t1637) * MDP(12) + (t1485 * t1642 + t1486 * t1640 - t1487 * t1638) * MDP(13)) * t1676; t1610 * MDP(4) + (-t1630 - t1631 - t1632) * MDP(5) + ((t1447 - t1663) * t1679 + (t1445 - t1661) * t1682 + (t1443 - t1659) * t1685) * MDP(6) + (t1479 * t1570 + t1480 * t1572 + t1481 * t1574) * MDP(12) + (t1479 * t1734 + t1480 * t1733 + t1481 * t1732) * MDP(13) + ((t1677 + t1680 + t1683) * MDP(11) * t1584 / t1583 + ((-t1557 * t1683 - t1558 * t1680 - t1559 * t1677) * MDP(12) + (t1678 + t1681 + t1684) * MDP(13)) * t1585) * t1560 + ((-t1643 - t1644 - t1645) * MDP(7) + (-t1488 * t1527 * t1693 - t1489 * t1528 * t1689 - t1490 * t1529 * t1686) * MDP(8) - t1610 * MDP(9) + (t1455 + t1456 + t1457) * MDP(10) + (qJ(2,1) * t1630 + qJ(2,2) * t1631 + qJ(2,3) * t1632 + t1440 * t1685 + t1441 * t1682 + t1442 * t1679) * MDP(12) + (-t1440 - t1441 - t1442 - t1659 - t1661 - t1663) * MDP(13)) * t1584;];
taucX  = t1;
