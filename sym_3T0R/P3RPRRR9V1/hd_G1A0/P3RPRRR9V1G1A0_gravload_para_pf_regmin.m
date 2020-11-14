% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x15]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR9V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:50
% EndTime: 2020-08-06 18:47:51
% DurationCPUTime: 0.82s
% Computational Cost: add. (729->140), mult. (915->245), div. (117->7), fcn. (954->29), ass. (0->104)
t1590 = legFrame(3,3);
t1573 = sin(t1590);
t1576 = cos(t1590);
t1557 = t1573 * g(1) - t1576 * g(2);
t1560 = t1576 * g(1) + t1573 * g(2);
t1593 = sin(qJ(1,3));
t1596 = cos(qJ(1,3));
t1530 = t1557 * t1593 - t1560 * t1596;
t1591 = legFrame(2,3);
t1574 = sin(t1591);
t1577 = cos(t1591);
t1558 = t1574 * g(1) - t1577 * g(2);
t1561 = t1577 * g(1) + t1574 * g(2);
t1594 = sin(qJ(1,2));
t1597 = cos(qJ(1,2));
t1532 = t1558 * t1594 - t1561 * t1597;
t1592 = legFrame(1,3);
t1575 = sin(t1592);
t1578 = cos(t1592);
t1559 = t1575 * g(1) - t1578 * g(2);
t1562 = t1578 * g(1) + t1575 * g(2);
t1595 = sin(qJ(1,1));
t1598 = cos(qJ(1,1));
t1534 = t1559 * t1595 - t1562 * t1598;
t1587 = pkin(7) + qJ(3,1);
t1572 = cos(t1587);
t1566 = 0.1e1 / t1572;
t1569 = sin(t1587);
t1631 = -pkin(5) - pkin(6);
t1584 = -qJ(2,1) + t1631;
t1581 = 0.1e1 / t1584;
t1616 = t1581 * t1569;
t1607 = t1566 * t1616;
t1586 = pkin(7) + qJ(3,2);
t1571 = cos(t1586);
t1565 = 0.1e1 / t1571;
t1568 = sin(t1586);
t1583 = -qJ(2,2) + t1631;
t1580 = 0.1e1 / t1583;
t1618 = t1580 * t1568;
t1608 = t1565 * t1618;
t1585 = pkin(7) + qJ(3,3);
t1570 = cos(t1585);
t1564 = 0.1e1 / t1570;
t1567 = sin(t1585);
t1582 = -qJ(2,3) + t1631;
t1579 = 0.1e1 / t1582;
t1620 = t1579 * t1567;
t1609 = t1564 * t1620;
t1635 = t1530 * t1609 + t1532 * t1608 + t1534 * t1607;
t1556 = t1575 * t1598 + t1578 * t1595;
t1622 = t1556 * t1581;
t1555 = t1574 * t1597 + t1577 * t1594;
t1623 = t1555 * t1580;
t1554 = t1573 * t1596 + t1576 * t1593;
t1624 = t1554 * t1579;
t1634 = t1530 * t1624 + t1532 * t1623 + t1534 * t1622;
t1604 = t1575 * t1595 - t1578 * t1598;
t1625 = t1604 * t1581;
t1605 = t1574 * t1594 - t1577 * t1597;
t1626 = t1605 * t1580;
t1606 = t1573 * t1593 - t1576 * t1596;
t1627 = t1606 * t1579;
t1633 = t1530 * t1627 + t1532 * t1626 + t1534 * t1625;
t1632 = 2 * pkin(1);
t1630 = pkin(3) * t1570;
t1629 = pkin(3) * t1571;
t1628 = pkin(3) * t1572;
t1527 = t1554 * g(1) + t1606 * g(2);
t1621 = t1579 * t1527;
t1528 = t1555 * g(1) + t1605 * g(2);
t1619 = t1580 * t1528;
t1529 = t1556 * g(1) + t1604 * g(2);
t1617 = t1581 * t1529;
t1615 = t1570 * t1621;
t1614 = t1571 * t1619;
t1613 = t1572 * t1617;
t1612 = t1527 * t1620;
t1611 = t1528 * t1618;
t1610 = t1529 * t1616;
t1531 = t1557 * t1596 + t1560 * t1593;
t1533 = t1558 * t1597 + t1561 * t1594;
t1535 = t1559 * t1598 + t1562 * t1595;
t1603 = -t1531 * t1627 - t1533 * t1626 - t1535 * t1625;
t1602 = t1531 * t1624 + t1533 * t1623 + t1535 * t1622;
t1601 = t1531 * t1609 + t1533 * t1608 + t1535 * t1607;
t1600 = 0.1e1 / pkin(3);
t1599 = 2 * pkin(7);
t1589 = cos(pkin(7));
t1588 = sin(pkin(7));
t1563 = t1589 * pkin(2) + pkin(1);
t1550 = t1563 * t1598 - t1595 * t1584;
t1549 = t1563 * t1597 - t1594 * t1583;
t1548 = t1563 * t1596 - t1593 * t1582;
t1547 = t1595 * t1563 + t1598 * t1584;
t1546 = t1594 * t1563 + t1597 * t1583;
t1545 = t1593 * t1563 + t1596 * t1582;
t1544 = (g(1) * t1598 + g(2) * t1595) * t1578 - (g(1) * t1595 - g(2) * t1598) * t1575;
t1543 = (t1597 * g(1) + t1594 * g(2)) * t1577 - (t1594 * g(1) - t1597 * g(2)) * t1574;
t1542 = (t1596 * g(1) + t1593 * g(2)) * t1576 - (t1593 * g(1) - t1596 * g(2)) * t1573;
t1526 = -t1562 * (-t1595 * pkin(1) + t1598 * qJ(2,1)) + t1559 * (t1598 * pkin(1) + t1595 * qJ(2,1));
t1525 = -t1561 * (-t1594 * pkin(1) + t1597 * qJ(2,2)) + t1558 * (t1597 * pkin(1) + t1594 * qJ(2,2));
t1524 = -t1560 * (-t1593 * pkin(1) + t1596 * qJ(2,3)) + t1557 * (t1596 * pkin(1) + t1593 * qJ(2,3));
t1 = [0, -t1603, -t1633, -t1603 * t1589, t1603 * t1588, t1633, -(-t1604 * t1526 - (-t1547 * t1575 + t1550 * t1578 - t1604 * t1628) * t1535) * t1581 - (-t1605 * t1525 - (-t1546 * t1574 + t1549 * t1577 - t1605 * t1629) * t1533) * t1580 - (-t1606 * t1524 - (-t1545 * t1573 + t1548 * t1576 - t1606 * t1630) * t1531) * t1579, 0, 0, 0, 0, 0, t1604 * t1613 + t1605 * t1614 + t1606 * t1615, -t1604 * t1610 - t1605 * t1611 - t1606 * t1612, -g(1); 0, -t1602, t1634, -t1602 * t1589, t1602 * t1588, -t1634, -(t1556 * t1526 - (t1547 * t1578 + t1550 * t1575 + t1556 * t1628) * t1535) * t1581 - (t1555 * t1525 - (t1546 * t1577 + t1549 * t1574 + t1555 * t1629) * t1533) * t1580 - (t1554 * t1524 - (t1545 * t1576 + t1548 * t1573 + t1554 * t1630) * t1531) * t1579, 0, 0, 0, 0, 0, -t1554 * t1615 - t1555 * t1614 - t1556 * t1613, t1554 * t1612 + t1555 * t1611 + t1556 * t1610, -g(2); 0, -t1601, t1635, -t1601 * t1589, t1601 * t1588, -t1635, -t1526 * t1607 + t1581 * (pkin(3) * sin(0.2e1 * t1587) + t1569 * t1632 + (sin((t1599 + qJ(3,1))) + sin(qJ(3,1))) * pkin(2)) * t1566 * t1535 / 0.2e1 - t1525 * t1608 + t1580 * (pkin(3) * sin(0.2e1 * t1586) + t1568 * t1632 + (sin((t1599 + qJ(3,2))) + sin(qJ(3,2))) * pkin(2)) * t1565 * t1533 / 0.2e1 - t1524 * t1609 + t1579 * (pkin(3) * sin(0.2e1 * t1585) + t1567 * t1632 + (sin((t1599 + qJ(3,3))) + sin(qJ(3,3))) * pkin(2)) * t1564 * t1531 / 0.2e1, 0, 0, 0, 0, 0, -t1612 - t1611 - t1610 + (t1566 * (-g(3) * t1572 + t1544 * t1569) + t1565 * (-g(3) * t1571 + t1543 * t1568) + t1564 * (-g(3) * t1570 + t1542 * t1567)) * t1600, t1567 ^ 2 * t1564 * t1621 + t1568 ^ 2 * t1565 * t1619 + t1569 ^ 2 * t1566 * t1617 + (t1566 * (g(3) * t1569 + t1544 * t1572) + t1565 * (g(3) * t1568 + t1543 * t1571) + t1564 * (g(3) * t1567 + t1542 * t1570)) * t1600, -g(3);];
tau_reg  = t1;
