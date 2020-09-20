% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:09:48
% EndTime: 2020-08-07 11:09:50
% DurationCPUTime: 2.27s
% Computational Cost: add. (1638->245), mult. (3336->459), div. (64->9), fcn. (2810->30), ass. (0->209)
t1595 = sin(qJ(2,4));
t1597 = cos(qJ(2,4));
t1590 = legFrame(4,3);
t1571 = sin(t1590);
t1575 = cos(t1590);
t1546 = -t1571 * g(1) + t1575 * g(2);
t1550 = t1575 * g(1) + t1571 * g(2);
t1586 = sin(pkin(8));
t1588 = cos(pkin(8));
t1644 = t1546 * t1586 + t1550 * t1588;
t1589 = cos(pkin(4));
t1645 = t1546 * t1588 - t1550 * t1586;
t1587 = sin(pkin(4));
t1694 = g(3) * t1587;
t1702 = t1645 * t1589 + t1694;
t1629 = t1702 * t1595 + t1644 * t1597;
t1599 = sin(qJ(2,3));
t1605 = cos(qJ(2,3));
t1591 = legFrame(3,3);
t1572 = sin(t1591);
t1576 = cos(t1591);
t1547 = -t1572 * g(1) + t1576 * g(2);
t1551 = t1576 * g(1) + t1572 * g(2);
t1642 = t1547 * t1586 + t1551 * t1588;
t1643 = t1547 * t1588 - t1551 * t1586;
t1703 = t1589 * t1643 + t1694;
t1628 = t1703 * t1599 + t1642 * t1605;
t1601 = sin(qJ(2,2));
t1607 = cos(qJ(2,2));
t1592 = legFrame(2,3);
t1573 = sin(t1592);
t1577 = cos(t1592);
t1548 = -t1573 * g(1) + t1577 * g(2);
t1552 = t1577 * g(1) + t1573 * g(2);
t1640 = t1548 * t1586 + t1552 * t1588;
t1641 = t1548 * t1588 - t1552 * t1586;
t1704 = t1589 * t1641 + t1694;
t1627 = t1704 * t1601 + t1640 * t1607;
t1603 = sin(qJ(2,1));
t1609 = cos(qJ(2,1));
t1593 = legFrame(1,3);
t1574 = sin(t1593);
t1578 = cos(t1593);
t1549 = -t1574 * g(1) + t1578 * g(2);
t1553 = t1578 * g(1) + t1574 * g(2);
t1638 = t1549 * t1586 + t1553 * t1588;
t1639 = t1549 * t1588 - t1553 * t1586;
t1705 = t1589 * t1639 + t1694;
t1626 = t1705 * t1603 + t1638 * t1609;
t1693 = g(3) * t1589;
t1612 = pkin(7) + pkin(6);
t1656 = t1603 * t1612;
t1561 = pkin(2) * t1609 + t1656;
t1570 = t1612 * t1609;
t1558 = pkin(2) * t1603 - t1570;
t1602 = sin(qJ(3,1));
t1671 = t1587 * t1602;
t1634 = pkin(3) * t1671 - t1558 * t1589;
t1709 = t1561 * t1588 + t1586 * t1634;
t1657 = t1601 * t1612;
t1560 = pkin(2) * t1607 + t1657;
t1569 = t1612 * t1607;
t1557 = pkin(2) * t1601 - t1569;
t1600 = sin(qJ(3,2));
t1672 = t1587 * t1600;
t1635 = pkin(3) * t1672 - t1557 * t1589;
t1708 = t1560 * t1588 + t1586 * t1635;
t1658 = t1599 * t1612;
t1559 = pkin(2) * t1605 + t1658;
t1568 = t1612 * t1605;
t1556 = pkin(2) * t1599 - t1568;
t1598 = sin(qJ(3,3));
t1673 = t1587 * t1598;
t1636 = pkin(3) * t1673 - t1556 * t1589;
t1707 = t1559 * t1588 + t1586 * t1636;
t1659 = t1595 * t1612;
t1555 = pkin(2) * t1597 + t1659;
t1564 = t1612 * t1597;
t1554 = pkin(2) * t1595 - t1564;
t1594 = sin(qJ(3,4));
t1675 = t1587 * t1594;
t1637 = pkin(3) * t1675 - t1554 * t1589;
t1706 = t1555 * t1588 + t1586 * t1637;
t1701 = m(3) / pkin(3);
t1596 = cos(qJ(3,4));
t1581 = t1596 ^ 2;
t1700 = pkin(3) * t1581;
t1604 = cos(qJ(3,3));
t1583 = t1604 ^ 2;
t1699 = pkin(3) * t1583;
t1606 = cos(qJ(3,2));
t1584 = t1606 ^ 2;
t1698 = pkin(3) * t1584;
t1608 = cos(qJ(3,1));
t1585 = t1608 ^ 2;
t1697 = pkin(3) * t1585;
t1696 = pkin(3) * t1587;
t1582 = m(1) + m(2) + m(3);
t1695 = g(3) * t1582;
t1692 = rSges(3,2) * t1587;
t1563 = t1596 * pkin(3) + pkin(2);
t1531 = t1595 * t1563 - t1564;
t1655 = rSges(3,2) * t1693;
t1667 = t1589 * t1594;
t1674 = t1587 * t1596;
t1691 = (((t1645 * t1587 - t1693) * rSges(3,1) + t1629 * rSges(3,2)) * t1596 + t1594 * (t1629 * rSges(3,1) - t1645 * t1692 + t1655)) / (t1531 * t1674 + t1563 * t1667);
t1565 = t1604 * pkin(3) + pkin(2);
t1535 = t1599 * t1565 - t1568;
t1665 = t1589 * t1598;
t1670 = t1587 * t1604;
t1690 = (((t1587 * t1643 - t1693) * rSges(3,1) + t1628 * rSges(3,2)) * t1604 + t1598 * (rSges(3,1) * t1628 - t1643 * t1692 + t1655)) / (t1535 * t1670 + t1565 * t1665);
t1566 = t1606 * pkin(3) + pkin(2);
t1536 = t1601 * t1566 - t1569;
t1663 = t1589 * t1600;
t1669 = t1587 * t1606;
t1689 = (((t1587 * t1641 - t1693) * rSges(3,1) + t1627 * rSges(3,2)) * t1606 + t1600 * (rSges(3,1) * t1627 - t1641 * t1692 + t1655)) / (t1536 * t1669 + t1566 * t1663);
t1567 = t1608 * pkin(3) + pkin(2);
t1537 = t1603 * t1567 - t1570;
t1661 = t1589 * t1602;
t1668 = t1587 * t1608;
t1688 = (((t1587 * t1639 - t1693) * rSges(3,1) + t1626 * rSges(3,2)) * t1608 + t1602 * (rSges(3,1) * t1626 - t1639 * t1692 + t1655)) / (t1537 * t1668 + t1567 * t1661);
t1494 = 0.1e1 / (t1595 * t1581 * t1696 + (pkin(3) * t1667 + t1554 * t1587) * t1596 + pkin(2) * t1667);
t1562 = (-pkin(6) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t1654 = m(2) * rSges(2,1) + pkin(2) * m(3);
t1687 = (t1629 * t1562 + (t1644 * t1595 - t1597 * t1702) * ((rSges(3,1) * t1596 - rSges(3,2) * t1594) * m(3) + t1654)) * t1494;
t1495 = 0.1e1 / (t1599 * t1583 * t1696 + (pkin(3) * t1665 + t1556 * t1587) * t1604 + pkin(2) * t1665);
t1686 = (t1628 * t1562 + (t1642 * t1599 - t1605 * t1703) * ((rSges(3,1) * t1604 - rSges(3,2) * t1598) * m(3) + t1654)) * t1495;
t1496 = 0.1e1 / (t1601 * t1584 * t1696 + (pkin(3) * t1663 + t1557 * t1587) * t1606 + pkin(2) * t1663);
t1685 = (t1627 * t1562 + (t1640 * t1601 - t1607 * t1704) * ((rSges(3,1) * t1606 - rSges(3,2) * t1600) * m(3) + t1654)) * t1496;
t1497 = 0.1e1 / (t1603 * t1585 * t1696 + (pkin(3) * t1661 + t1558 * t1587) * t1608 + pkin(2) * t1661);
t1684 = (t1626 * t1562 + (t1638 * t1603 - t1609 * t1705) * ((rSges(3,1) * t1608 - rSges(3,2) * t1602) * m(3) + t1654)) * t1497;
t1683 = (t1563 * t1597 + t1659) * t1589;
t1682 = (t1565 * t1605 + t1658) * t1589;
t1681 = (t1566 * t1607 + t1657) * t1589;
t1680 = (t1567 * t1609 + t1656) * t1589;
t1666 = t1589 * t1595;
t1664 = t1589 * t1599;
t1662 = t1589 * t1601;
t1660 = t1589 * t1603;
t1653 = pkin(2) * t1675;
t1652 = pkin(2) * t1673;
t1651 = pkin(2) * t1672;
t1650 = pkin(2) * t1671;
t1624 = koppelP(1,1);
t1623 = koppelP(2,1);
t1622 = koppelP(3,1);
t1621 = koppelP(4,1);
t1620 = koppelP(1,2);
t1619 = koppelP(2,2);
t1618 = koppelP(3,2);
t1617 = koppelP(4,2);
t1616 = rSges(4,1);
t1615 = rSges(4,2);
t1614 = xP(4);
t1580 = cos(t1614);
t1579 = sin(t1614);
t1545 = -t1579 * t1620 + t1580 * t1624;
t1544 = -t1579 * t1619 + t1580 * t1623;
t1543 = -t1579 * t1618 + t1580 * t1622;
t1542 = -t1579 * t1617 + t1580 * t1621;
t1541 = -t1579 * t1624 - t1580 * t1620;
t1540 = -t1579 * t1623 - t1580 * t1619;
t1539 = -t1579 * t1622 - t1580 * t1618;
t1538 = -t1579 * t1621 - t1580 * t1617;
t1529 = t1586 * t1609 + t1588 * t1660;
t1528 = t1586 * t1607 + t1588 * t1662;
t1527 = t1586 * t1605 + t1588 * t1664;
t1526 = t1586 * t1660 - t1588 * t1609;
t1525 = t1586 * t1662 - t1588 * t1607;
t1524 = t1586 * t1664 - t1588 * t1605;
t1520 = t1586 * t1597 + t1588 * t1666;
t1519 = t1586 * t1666 - t1588 * t1597;
t1517 = t1588 * t1574 + t1578 * t1586;
t1516 = t1588 * t1573 + t1577 * t1586;
t1515 = t1588 * t1572 + t1576 * t1586;
t1514 = t1588 * t1571 + t1575 * t1586;
t1513 = -t1586 * t1574 + t1578 * t1588;
t1512 = -t1586 * t1573 + t1577 * t1588;
t1511 = -t1586 * t1572 + t1576 * t1588;
t1510 = -t1586 * t1571 + t1575 * t1588;
t1501 = t1586 * t1561 - t1588 * t1634;
t1500 = t1586 * t1560 - t1588 * t1635;
t1499 = t1586 * t1559 - t1588 * t1636;
t1498 = t1586 * t1555 - t1588 * t1637;
t1493 = -t1537 * t1513 - t1517 * t1680;
t1492 = -t1536 * t1512 - t1516 * t1681;
t1491 = -t1535 * t1511 - t1515 * t1682;
t1490 = -t1513 * t1680 + t1537 * t1517;
t1489 = -t1512 * t1681 + t1536 * t1516;
t1488 = -t1511 * t1682 + t1535 * t1515;
t1487 = -t1531 * t1510 - t1514 * t1683;
t1486 = -t1510 * t1683 + t1531 * t1514;
t1485 = -t1517 * t1668 - (-t1609 * t1513 + t1517 * t1660) * t1602;
t1484 = -t1516 * t1669 - (-t1607 * t1512 + t1516 * t1662) * t1600;
t1483 = -t1515 * t1670 - (-t1605 * t1511 + t1515 * t1664) * t1598;
t1482 = -t1513 * t1668 - (t1513 * t1660 + t1609 * t1517) * t1602;
t1481 = -t1512 * t1669 - (t1512 * t1662 + t1607 * t1516) * t1600;
t1480 = -t1511 * t1670 - (t1511 * t1664 + t1605 * t1515) * t1598;
t1479 = -t1514 * t1674 - (-t1597 * t1510 + t1514 * t1666) * t1594;
t1478 = -t1510 * t1674 - (t1510 * t1666 + t1597 * t1514) * t1594;
t1477 = (-t1574 * t1526 + t1529 * t1578) * t1697 + (t1501 * t1578 + t1709 * t1574) * t1608 - t1513 * t1650;
t1476 = (-t1573 * t1525 + t1528 * t1577) * t1698 + (t1500 * t1577 + t1708 * t1573) * t1606 - t1512 * t1651;
t1475 = (-t1572 * t1524 + t1527 * t1576) * t1699 + (t1499 * t1576 + t1707 * t1572) * t1604 - t1511 * t1652;
t1474 = -(t1526 * t1578 + t1574 * t1529) * t1697 + (-t1574 * t1501 + t1709 * t1578) * t1608 + t1517 * t1650;
t1473 = -(t1525 * t1577 + t1573 * t1528) * t1698 + (-t1573 * t1500 + t1708 * t1577) * t1606 + t1516 * t1651;
t1472 = -(t1524 * t1576 + t1572 * t1527) * t1699 + (-t1572 * t1499 + t1707 * t1576) * t1604 + t1515 * t1652;
t1471 = (-t1571 * t1519 + t1520 * t1575) * t1700 + (t1498 * t1575 + t1706 * t1571) * t1596 - t1510 * t1653;
t1470 = -(t1519 * t1575 + t1571 * t1520) * t1700 + (-t1571 * t1498 + t1706 * t1575) * t1596 + t1514 * t1653;
t1 = [t1478 * t1687 + t1480 * t1686 + t1481 * t1685 + t1482 * t1684 - m(4) * g(1) + (-t1470 * t1494 - t1472 * t1495 - t1473 * t1496 - t1474 * t1497) * t1695 + (t1486 * t1691 + t1488 * t1690 + t1489 * t1689 + t1490 * t1688) * t1701; t1479 * t1687 + t1483 * t1686 + t1484 * t1685 + t1485 * t1684 - m(4) * g(2) + (-t1471 * t1494 - t1475 * t1495 - t1476 * t1496 - t1477 * t1497) * t1695 + (t1487 * t1691 + t1491 * t1690 + t1492 * t1689 + t1493 * t1688) * t1701; (-m(4) - 0.4e1 * t1582) * g(3); m(4) * ((g(1) * t1616 + g(2) * t1615) * t1579 + (g(1) * t1615 - g(2) * t1616) * t1580) + (t1482 * t1541 + t1485 * t1545) * t1684 + (t1481 * t1540 + t1484 * t1544) * t1685 + (t1480 * t1539 + t1483 * t1543) * t1686 + (t1478 * t1538 + t1479 * t1542) * t1687 + ((-t1474 * t1541 - t1477 * t1545) * t1497 + (-t1473 * t1540 - t1476 * t1544) * t1496 + (-t1472 * t1539 - t1475 * t1543) * t1495 + (-t1470 * t1538 - t1471 * t1542) * t1494) * t1695 + ((t1490 * t1541 + t1493 * t1545) * t1688 + (t1489 * t1540 + t1492 * t1544) * t1689 + (t1488 * t1539 + t1491 * t1543) * t1690 + (t1486 * t1538 + t1487 * t1542) * t1691) * t1701;];
taugX  = t1;
