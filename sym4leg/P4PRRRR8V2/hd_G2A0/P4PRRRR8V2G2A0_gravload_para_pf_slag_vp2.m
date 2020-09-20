% Calculate Gravitation load for parallel robot
% P4PRRRR8V2G2A0
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
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [4x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:17:18
% EndTime: 2020-08-07 11:17:21
% DurationCPUTime: 2.80s
% Computational Cost: add. (1851->256), mult. (3705->488), div. (80->5), fcn. (3150->30), ass. (0->202)
t1569 = sin(qJ(2,4));
t1571 = cos(qJ(2,4));
t1572 = legFrame(4,2);
t1549 = sin(t1572);
t1553 = cos(t1572);
t1534 = t1553 * g(1) - t1549 * g(2);
t1564 = sin(pkin(8));
t1566 = cos(pkin(8));
t1692 = g(3) * t1566;
t1625 = t1534 * t1564 + t1692;
t1530 = t1549 * g(1) + t1553 * g(2);
t1565 = sin(pkin(4));
t1567 = cos(pkin(4));
t1548 = g(3) * t1564;
t1624 = t1534 * t1566 - t1548;
t1713 = t1530 * t1565 + t1624 * t1567;
t1605 = t1713 * t1569 + t1625 * t1571;
t1577 = sin(qJ(2,3));
t1583 = cos(qJ(2,3));
t1573 = legFrame(3,2);
t1550 = sin(t1573);
t1554 = cos(t1573);
t1535 = t1554 * g(1) - t1550 * g(2);
t1623 = t1535 * t1564 + t1692;
t1531 = t1550 * g(1) + t1554 * g(2);
t1622 = t1535 * t1566 - t1548;
t1714 = t1531 * t1565 + t1622 * t1567;
t1604 = t1714 * t1577 + t1623 * t1583;
t1579 = sin(qJ(2,2));
t1585 = cos(qJ(2,2));
t1574 = legFrame(2,2);
t1551 = sin(t1574);
t1555 = cos(t1574);
t1536 = t1555 * g(1) - t1551 * g(2);
t1621 = t1536 * t1564 + t1692;
t1532 = t1551 * g(1) + t1555 * g(2);
t1620 = t1536 * t1566 - t1548;
t1715 = t1532 * t1565 + t1620 * t1567;
t1603 = t1715 * t1579 + t1621 * t1585;
t1581 = sin(qJ(2,1));
t1587 = cos(qJ(2,1));
t1575 = legFrame(1,2);
t1552 = sin(t1575);
t1556 = cos(t1575);
t1537 = t1556 * g(1) - t1552 * g(2);
t1619 = t1537 * t1564 + t1692;
t1533 = t1552 * g(1) + t1556 * g(2);
t1618 = t1537 * t1566 - t1548;
t1716 = t1533 * t1565 + t1618 * t1567;
t1602 = t1716 * t1581 + t1619 * t1587;
t1568 = sin(qJ(3,4));
t1704 = pkin(2) * t1568;
t1576 = sin(qJ(3,3));
t1703 = pkin(2) * t1576;
t1578 = sin(qJ(3,2));
t1702 = pkin(2) * t1578;
t1580 = sin(qJ(3,1));
t1701 = pkin(2) * t1580;
t1570 = cos(qJ(3,4));
t1700 = pkin(3) * t1570 ^ 2;
t1582 = cos(qJ(3,3));
t1699 = pkin(3) * t1582 ^ 2;
t1584 = cos(qJ(3,2));
t1698 = pkin(3) * t1584 ^ 2;
t1586 = cos(qJ(3,1));
t1697 = pkin(3) * t1586 ^ 2;
t1696 = pkin(3) * t1570;
t1695 = pkin(3) * t1582;
t1694 = pkin(3) * t1584;
t1693 = pkin(3) * t1586;
t1691 = m(3) * pkin(2) + mrSges(2,1);
t1588 = pkin(7) + pkin(6);
t1538 = pkin(2) * t1569 - t1588 * t1571;
t1645 = t1567 * t1568;
t1502 = pkin(3) * t1645 + t1538 * t1565;
t1657 = t1565 * t1569;
t1490 = 0.1e1 / (pkin(2) * t1645 + t1502 * t1570 + t1657 * t1700);
t1659 = t1565 * t1566;
t1662 = mrSges(3,2) * t1548 * t1565;
t1676 = t1530 * t1567;
t1690 = (((t1624 * t1565 - t1676) * mrSges(3,1) + t1605 * mrSges(3,2)) * t1570 + (t1662 + (-t1534 * t1659 + t1676) * mrSges(3,2) + t1605 * mrSges(3,1)) * t1568) * t1490;
t1540 = pkin(2) * t1577 - t1588 * t1583;
t1642 = t1567 * t1576;
t1503 = pkin(3) * t1642 + t1540 * t1565;
t1655 = t1565 * t1577;
t1491 = 0.1e1 / (pkin(2) * t1642 + t1503 * t1582 + t1655 * t1699);
t1673 = t1531 * t1567;
t1689 = (((t1622 * t1565 - t1673) * mrSges(3,1) + t1604 * mrSges(3,2)) * t1582 + (t1662 + (-t1535 * t1659 + t1673) * mrSges(3,2) + t1604 * mrSges(3,1)) * t1576) * t1491;
t1541 = pkin(2) * t1579 - t1588 * t1585;
t1640 = t1567 * t1578;
t1504 = pkin(3) * t1640 + t1541 * t1565;
t1653 = t1565 * t1579;
t1492 = 0.1e1 / (pkin(2) * t1640 + t1504 * t1584 + t1653 * t1698);
t1670 = t1532 * t1567;
t1688 = (((t1620 * t1565 - t1670) * mrSges(3,1) + t1603 * mrSges(3,2)) * t1584 + (t1662 + (-t1536 * t1659 + t1670) * mrSges(3,2) + t1603 * mrSges(3,1)) * t1578) * t1492;
t1542 = pkin(2) * t1581 - t1588 * t1587;
t1638 = t1567 * t1580;
t1505 = pkin(3) * t1638 + t1542 * t1565;
t1651 = t1565 * t1581;
t1493 = 0.1e1 / (pkin(2) * t1638 + t1505 * t1586 + t1651 * t1697);
t1667 = t1533 * t1567;
t1687 = (((t1618 * t1565 - t1667) * mrSges(3,1) + t1602 * mrSges(3,2)) * t1586 + (t1662 + (-t1537 * t1659 + t1667) * mrSges(3,2) + t1602 * mrSges(3,1)) * t1580) * t1493;
t1547 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t1474 = -t1605 * t1547 + (t1625 * t1569 - t1571 * t1713) * (t1570 * mrSges(3,1) - mrSges(3,2) * t1568 + t1691);
t1644 = t1567 * t1569;
t1510 = t1564 * t1644 - t1566 * t1571;
t1661 = t1564 * t1565;
t1686 = t1474 * (t1568 * t1510 + t1570 * t1661);
t1475 = -t1604 * t1547 + (t1623 * t1577 - t1583 * t1714) * (t1582 * mrSges(3,1) - mrSges(3,2) * t1576 + t1691);
t1641 = t1567 * t1577;
t1512 = t1564 * t1641 - t1566 * t1583;
t1685 = t1475 * (t1576 * t1512 + t1582 * t1661);
t1476 = -t1603 * t1547 + (t1621 * t1579 - t1585 * t1715) * (t1584 * mrSges(3,1) - mrSges(3,2) * t1578 + t1691);
t1639 = t1567 * t1579;
t1513 = t1564 * t1639 - t1566 * t1585;
t1684 = t1476 * (t1578 * t1513 + t1584 * t1661);
t1477 = -t1602 * t1547 + (t1619 * t1581 - t1587 * t1716) * (t1586 * mrSges(3,1) - mrSges(3,2) * t1580 + t1691);
t1637 = t1567 * t1581;
t1514 = t1564 * t1637 - t1566 * t1587;
t1683 = t1477 * (t1580 * t1514 + t1586 * t1661);
t1682 = t1490 * t1530;
t1681 = t1491 * t1531;
t1680 = t1492 * t1532;
t1679 = t1493 * t1533;
t1560 = m(1) + m(2) + m(3);
t1678 = t1530 * t1560;
t1675 = t1531 * t1560;
t1672 = t1532 * t1560;
t1669 = t1533 * t1560;
t1660 = t1564 * t1567;
t1658 = t1565 * t1568;
t1656 = t1565 * t1576;
t1654 = t1565 * t1578;
t1652 = t1565 * t1580;
t1650 = t1566 * t1567;
t1649 = t1566 * t1570;
t1648 = t1566 * t1582;
t1647 = t1566 * t1584;
t1646 = t1566 * t1586;
t1643 = t1567 * t1571;
t1636 = t1567 * t1583;
t1635 = t1567 * t1585;
t1634 = t1567 * t1587;
t1539 = pkin(2) * t1571 + t1569 * t1588;
t1633 = ((t1564 * t1643 + t1566 * t1569) * t1696 + t1539 * t1660 + t1538 * t1566) * t1690;
t1543 = pkin(2) * t1583 + t1577 * t1588;
t1632 = ((t1564 * t1636 + t1566 * t1577) * t1695 + t1543 * t1660 + t1540 * t1566) * t1689;
t1544 = pkin(2) * t1585 + t1579 * t1588;
t1631 = ((t1564 * t1635 + t1566 * t1579) * t1694 + t1544 * t1660 + t1541 * t1566) * t1688;
t1545 = pkin(2) * t1587 + t1581 * t1588;
t1630 = ((t1564 * t1634 + t1566 * t1581) * t1693 + t1545 * t1660 + t1542 * t1566) * t1687;
t1629 = t1490 * t1686;
t1628 = t1491 * t1685;
t1627 = t1492 * t1684;
t1626 = t1493 * t1683;
t1590 = xP(4);
t1557 = sin(t1590);
t1558 = cos(t1590);
t1593 = koppelP(4,2);
t1597 = koppelP(4,1);
t1518 = -t1557 * t1597 - t1558 * t1593;
t1522 = -t1557 * t1593 + t1558 * t1597;
t1613 = -t1518 * t1553 + t1522 * t1549;
t1594 = koppelP(3,2);
t1598 = koppelP(3,1);
t1519 = -t1557 * t1598 - t1558 * t1594;
t1523 = -t1557 * t1594 + t1558 * t1598;
t1612 = -t1519 * t1554 + t1523 * t1550;
t1595 = koppelP(2,2);
t1599 = koppelP(2,1);
t1520 = -t1557 * t1599 - t1558 * t1595;
t1524 = -t1557 * t1595 + t1558 * t1599;
t1611 = -t1520 * t1555 + t1524 * t1551;
t1596 = koppelP(1,2);
t1600 = koppelP(1,1);
t1521 = -t1557 * t1600 - t1558 * t1596;
t1525 = -t1557 * t1596 + t1558 * t1600;
t1610 = -t1521 * t1556 + t1525 * t1552;
t1609 = pkin(3) * t1658 - t1538 * t1567;
t1608 = pkin(3) * t1656 - t1540 * t1567;
t1607 = pkin(3) * t1654 - t1541 * t1567;
t1606 = pkin(3) * t1652 - t1542 * t1567;
t1601 = 0.1e1 / pkin(3);
t1592 = mrSges(4,1);
t1591 = mrSges(4,2);
t1517 = t1564 * t1587 + t1566 * t1637;
t1516 = t1564 * t1585 + t1566 * t1639;
t1515 = t1564 * t1583 + t1566 * t1641;
t1511 = t1564 * t1571 + t1566 * t1644;
t1497 = -t1545 * t1564 + t1606 * t1566;
t1496 = -t1544 * t1564 + t1607 * t1566;
t1495 = -t1543 * t1564 + t1608 * t1566;
t1494 = -t1539 * t1564 + t1609 * t1566;
t1485 = -(t1517 * t1552 - t1556 * t1651) * t1697 + (t1497 * t1552 + t1556 * t1505) * t1586 + (t1552 * t1659 + t1567 * t1556) * t1701;
t1484 = -(t1516 * t1551 - t1555 * t1653) * t1698 + (t1496 * t1551 + t1555 * t1504) * t1584 + (t1551 * t1659 + t1567 * t1555) * t1702;
t1483 = -(t1515 * t1550 - t1554 * t1655) * t1699 + (t1495 * t1550 + t1554 * t1503) * t1582 + (t1550 * t1659 + t1567 * t1554) * t1703;
t1482 = (t1517 * t1556 + t1552 * t1651) * t1697 + (-t1497 * t1556 + t1552 * t1505) * t1586 + (t1567 * t1552 - t1556 * t1659) * t1701;
t1481 = (t1516 * t1555 + t1551 * t1653) * t1698 + (-t1496 * t1555 + t1551 * t1504) * t1584 + (t1567 * t1551 - t1555 * t1659) * t1702;
t1480 = (t1515 * t1554 + t1550 * t1655) * t1699 + (-t1495 * t1554 + t1550 * t1503) * t1582 + (t1567 * t1550 - t1554 * t1659) * t1703;
t1479 = -(t1511 * t1549 - t1553 * t1657) * t1700 + (t1494 * t1549 + t1553 * t1502) * t1570 + (t1549 * t1659 + t1567 * t1553) * t1704;
t1478 = (t1511 * t1553 + t1549 * t1657) * t1700 + (-t1494 * t1553 + t1549 * t1502) * t1570 + (t1567 * t1549 - t1553 * t1659) * t1704;
t1 = [-t1553 * t1629 - t1554 * t1628 - t1555 * t1627 - t1556 * t1626 - g(1) * m(4) + (-t1553 * t1633 - t1554 * t1632 - t1555 * t1631 - t1556 * t1630) * t1601 + (-t1478 * t1682 - t1480 * t1681 - t1481 * t1680 - t1482 * t1679) * t1560; t1549 * t1629 + t1550 * t1628 + t1551 * t1627 + t1552 * t1626 - g(2) * m(4) + (t1549 * t1633 + t1550 * t1632 + t1551 * t1631 + t1552 * t1630) * t1601 + (-t1479 * t1682 - t1483 * t1681 - t1484 * t1680 - t1485 * t1679) * t1560; -g(3) * m(4) + (-(-t1514 * t1697 + t1545 * t1646 + (pkin(2) * t1652 + t1606 * t1586) * t1564) * t1669 + (-t1580 * t1517 - t1565 * t1646) * t1477) * t1493 + (-(-t1513 * t1698 + t1544 * t1647 + (pkin(2) * t1654 + t1607 * t1584) * t1564) * t1672 + (-t1578 * t1516 - t1565 * t1647) * t1476) * t1492 + (-(-t1512 * t1699 + t1543 * t1648 + (pkin(2) * t1656 + t1608 * t1582) * t1564) * t1675 + (-t1576 * t1515 - t1565 * t1648) * t1475) * t1491 + (-(-t1510 * t1700 + t1539 * t1649 + (pkin(2) * t1658 + t1609 * t1570) * t1564) * t1678 + (-t1568 * t1511 - t1565 * t1649) * t1474) * t1490 + (((t1564 * t1581 - t1566 * t1634) * t1693 - t1545 * t1650 + t1542 * t1564) * t1687 + ((t1564 * t1579 - t1566 * t1635) * t1694 - t1544 * t1650 + t1541 * t1564) * t1688 + ((t1564 * t1577 - t1566 * t1636) * t1695 - t1543 * t1650 + t1540 * t1564) * t1689 + ((t1564 * t1569 - t1566 * t1643) * t1696 - t1539 * t1650 + t1538 * t1564) * t1690) * t1601; -(-g(1) * t1592 - g(2) * t1591) * t1557 + t1558 * (g(1) * t1591 - g(2) * t1592) + (-(t1482 * t1521 + t1485 * t1525) * t1669 + t1610 * t1683) * t1493 + (-(t1481 * t1520 + t1484 * t1524) * t1672 + t1611 * t1684) * t1492 + (-(t1480 * t1519 + t1483 * t1523) * t1675 + t1612 * t1685) * t1491 + (-(t1478 * t1518 + t1479 * t1522) * t1678 + t1613 * t1686) * t1490 + (t1610 * t1630 + t1611 * t1631 + t1612 * t1632 + t1613 * t1633) * t1601;];
taugX  = t1;
