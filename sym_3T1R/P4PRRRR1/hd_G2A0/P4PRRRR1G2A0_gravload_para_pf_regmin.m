% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4PRRRR1G2A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR1G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:59:23
% EndTime: 2020-08-07 10:59:24
% DurationCPUTime: 1.05s
% Computational Cost: add. (433->121), mult. (1052->268), div. (272->13), fcn. (1176->26), ass. (0->139)
t1573 = legFrame(4,2);
t1540 = sin(t1573);
t1544 = cos(t1573);
t1529 = g(1) * t1540 + g(2) * t1544;
t1570 = sin(qJ(2,4));
t1572 = cos(qJ(2,4));
t1506 = -g(3) * t1570 + t1529 * t1572;
t1554 = 0.1e1 / t1570;
t1571 = cos(qJ(3,4));
t1555 = 0.1e1 / t1571;
t1658 = t1554 * t1555;
t1670 = t1506 * t1658;
t1574 = legFrame(3,2);
t1541 = sin(t1574);
t1545 = cos(t1574);
t1530 = g(1) * t1541 + g(2) * t1545;
t1578 = sin(qJ(2,3));
t1584 = cos(qJ(2,3));
t1510 = -g(3) * t1578 + t1530 * t1584;
t1558 = 0.1e1 / t1578;
t1583 = cos(qJ(3,3));
t1563 = 0.1e1 / t1583;
t1656 = t1558 * t1563;
t1669 = t1510 * t1656;
t1575 = legFrame(2,2);
t1542 = sin(t1575);
t1546 = cos(t1575);
t1531 = g(1) * t1542 + g(2) * t1546;
t1580 = sin(qJ(2,2));
t1586 = cos(qJ(2,2));
t1511 = -g(3) * t1580 + t1531 * t1586;
t1560 = 0.1e1 / t1580;
t1585 = cos(qJ(3,2));
t1565 = 0.1e1 / t1585;
t1654 = t1560 * t1565;
t1668 = t1511 * t1654;
t1576 = legFrame(1,2);
t1543 = sin(t1576);
t1547 = cos(t1576);
t1532 = g(1) * t1543 + g(2) * t1547;
t1582 = sin(qJ(2,1));
t1588 = cos(qJ(2,1));
t1512 = -g(3) * t1582 + t1532 * t1588;
t1562 = 0.1e1 / t1582;
t1587 = cos(qJ(3,1));
t1567 = 0.1e1 / t1587;
t1652 = t1562 * t1567;
t1667 = t1512 * t1652;
t1536 = g(1) * t1547 - g(2) * t1543;
t1581 = sin(qJ(3,1));
t1635 = g(3) * t1588 + t1532 * t1582;
t1651 = t1562 * t1588;
t1639 = t1581 * t1651;
t1603 = (t1512 * t1639 + t1536 * t1587 + t1581 * t1635) * t1567;
t1535 = g(1) * t1546 - g(2) * t1542;
t1579 = sin(qJ(3,2));
t1636 = g(3) * t1586 + t1531 * t1580;
t1653 = t1560 * t1586;
t1640 = t1579 * t1653;
t1604 = (t1511 * t1640 + t1535 * t1585 + t1579 * t1636) * t1565;
t1534 = g(1) * t1545 - g(2) * t1541;
t1577 = sin(qJ(3,3));
t1637 = g(3) * t1584 + t1530 * t1578;
t1655 = t1558 * t1584;
t1641 = t1577 * t1655;
t1605 = (t1510 * t1641 + t1534 * t1583 + t1577 * t1637) * t1563;
t1533 = g(1) * t1544 - g(2) * t1540;
t1569 = sin(qJ(3,4));
t1638 = g(3) * t1572 + t1529 * t1570;
t1657 = t1554 * t1572;
t1642 = t1569 * t1657;
t1606 = (t1506 * t1642 + t1533 * t1571 + t1569 * t1638) * t1555;
t1650 = t1570 * t1571;
t1649 = t1578 * t1583;
t1648 = t1580 * t1585;
t1647 = t1582 * t1587;
t1646 = t1529 * t1658;
t1645 = t1530 * t1656;
t1644 = t1531 * t1654;
t1643 = t1532 * t1652;
t1556 = 0.1e1 / t1571 ^ 2;
t1634 = t1556 * t1642;
t1564 = 0.1e1 / t1583 ^ 2;
t1633 = t1564 * t1641;
t1566 = 0.1e1 / t1585 ^ 2;
t1632 = t1566 * t1640;
t1568 = 0.1e1 / t1587 ^ 2;
t1631 = t1568 * t1639;
t1589 = xP(4);
t1551 = sin(t1589);
t1552 = cos(t1589);
t1593 = koppelP(1,2);
t1597 = koppelP(1,1);
t1611 = -t1551 * t1593 + t1552 * t1597;
t1612 = t1551 * t1597 + t1552 * t1593;
t1493 = t1611 * t1543 + t1612 * t1547;
t1630 = t1493 * t1631;
t1590 = koppelP(4,2);
t1594 = koppelP(4,1);
t1617 = -t1551 * t1590 + t1552 * t1594;
t1618 = t1551 * t1594 + t1552 * t1590;
t1494 = t1617 * t1540 + t1618 * t1544;
t1629 = t1494 * t1634;
t1591 = koppelP(3,2);
t1595 = koppelP(3,1);
t1615 = -t1551 * t1591 + t1552 * t1595;
t1616 = t1551 * t1595 + t1552 * t1591;
t1495 = t1615 * t1541 + t1616 * t1545;
t1628 = t1495 * t1633;
t1592 = koppelP(2,2);
t1596 = koppelP(2,1);
t1613 = -t1551 * t1592 + t1552 * t1596;
t1614 = t1551 * t1596 + t1552 * t1592;
t1496 = t1613 * t1542 + t1614 * t1546;
t1627 = t1496 * t1632;
t1626 = t1540 * t1634;
t1625 = t1541 * t1633;
t1624 = t1542 * t1632;
t1623 = t1543 * t1631;
t1622 = t1544 * t1634;
t1621 = t1545 * t1633;
t1620 = t1546 * t1632;
t1619 = t1547 * t1631;
t1602 = t1506 * t1569 ^ 2 * t1556 * t1657 - (-t1533 * t1569 + t1571 * t1638) * t1555;
t1601 = t1510 * t1577 ^ 2 * t1564 * t1655 - (-t1534 * t1577 + t1583 * t1637) * t1563;
t1600 = t1511 * t1579 ^ 2 * t1566 * t1653 - (-t1535 * t1579 + t1585 * t1636) * t1565;
t1599 = t1512 * t1581 ^ 2 * t1568 * t1651 - (-t1536 * t1581 + t1587 * t1635) * t1567;
t1598 = 0.1e1 / pkin(2);
t1538 = g(1) * t1552 + g(2) * t1551;
t1537 = g(1) * t1551 - g(2) * t1552;
t1528 = t1543 * t1581 + t1547 * t1647;
t1527 = t1542 * t1579 + t1546 * t1648;
t1526 = t1541 * t1577 + t1545 * t1649;
t1525 = t1543 * t1647 - t1547 * t1581;
t1524 = t1542 * t1648 - t1546 * t1579;
t1523 = t1541 * t1649 - t1545 * t1577;
t1522 = t1540 * t1569 + t1544 * t1650;
t1521 = t1540 * t1650 - t1544 * t1569;
t1 = [-t1521 * t1646 - t1523 * t1645 - t1524 * t1644 - t1525 * t1643, 0, (-t1506 * t1622 - t1510 * t1621 - t1511 * t1620 - t1512 * t1619) * t1598, (t1619 * t1635 + t1620 * t1636 + t1621 * t1637 + t1622 * t1638) * t1598, 0, 0, 0, 0, 0, (-t1544 * t1606 - t1545 * t1605 - t1546 * t1604 - t1547 * t1603) * t1598, (t1602 * t1544 + t1601 * t1545 + t1600 * t1546 + t1599 * t1547) * t1598, 0, 0, 0, -t1537 * t1551 - t1538 * t1552; -t1522 * t1646 - t1526 * t1645 - t1527 * t1644 - t1528 * t1643, 0, (t1506 * t1626 + t1510 * t1625 + t1511 * t1624 + t1512 * t1623) * t1598, (-t1623 * t1635 - t1624 * t1636 - t1625 * t1637 - t1626 * t1638) * t1598, 0, 0, 0, 0, 0, (t1540 * t1606 + t1541 * t1605 + t1542 * t1604 + t1543 * t1603) * t1598, (-t1602 * t1540 - t1601 * t1541 - t1600 * t1542 - t1599 * t1543) * t1598, 0, 0, 0, t1537 * t1552 - t1538 * t1551; -t1529 * t1657 - t1530 * t1655 - t1531 * t1653 - t1532 * t1651, 0, (t1667 + t1668 + t1669 + t1670) * t1598, (-t1635 * t1652 - t1636 * t1654 - t1637 * t1656 - t1638 * t1658) * t1598, 0, 0, 0, 0, 0, (t1506 * t1554 + t1510 * t1558 + t1511 * t1560 + t1512 * t1562) * t1598, (-t1569 * t1670 - t1577 * t1669 - t1579 * t1668 - t1581 * t1667) * t1598, 0, 0, 0, -g(3); -(-t1525 * t1612 + t1528 * t1611) * t1643 - (-t1524 * t1614 + t1527 * t1613) * t1644 - (-t1523 * t1616 + t1526 * t1615) * t1645 - (-t1521 * t1618 + t1522 * t1617) * t1646, 0, (t1506 * t1629 + t1510 * t1628 + t1511 * t1627 + t1512 * t1630) * t1598, (-t1627 * t1636 - t1628 * t1637 - t1629 * t1638 - t1630 * t1635) * t1598, 0, 0, 0, 0, 0, (t1493 * t1603 + t1494 * t1606 + t1495 * t1605 + t1496 * t1604) * t1598, (-t1599 * t1493 - t1602 * t1494 - t1601 * t1495 - t1600 * t1496) * t1598, 0, t1537, t1538, 0;];
tau_reg  = t1;
