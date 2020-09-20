% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G3P3A0
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRR1G3P3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G3P3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:07:05
% EndTime: 2020-03-09 21:07:06
% DurationCPUTime: 0.46s
% Computational Cost: add. (876->93), mult. (485->167), div. (126->5), fcn. (558->18), ass. (0->78)
t1600 = legFrame(3,2);
t1591 = sin(t1600);
t1597 = pkin(7) + qJ(2,3);
t1588 = qJ(3,3) + t1597;
t1576 = sin(t1588);
t1579 = cos(t1588);
t1582 = sin(t1597);
t1585 = cos(t1597);
t1625 = 0.1e1 / (t1576 * t1585 - t1582 * t1579);
t1631 = t1591 * t1625;
t1601 = legFrame(2,2);
t1592 = sin(t1601);
t1598 = pkin(7) + qJ(2,2);
t1589 = qJ(3,2) + t1598;
t1577 = sin(t1589);
t1580 = cos(t1589);
t1583 = sin(t1598);
t1586 = cos(t1598);
t1624 = 0.1e1 / (t1577 * t1586 - t1583 * t1580);
t1630 = t1592 * t1624;
t1602 = legFrame(1,2);
t1593 = sin(t1602);
t1599 = pkin(7) + qJ(2,1);
t1590 = qJ(3,1) + t1599;
t1578 = sin(t1590);
t1581 = cos(t1590);
t1584 = sin(t1599);
t1587 = cos(t1599);
t1623 = 0.1e1 / (t1578 * t1587 - t1584 * t1581);
t1629 = t1593 * t1623;
t1594 = cos(t1600);
t1628 = t1594 * t1625;
t1595 = cos(t1601);
t1627 = t1595 * t1624;
t1596 = cos(t1602);
t1626 = t1596 * t1623;
t1622 = t1625 * (pkin(2) * t1582 + pkin(3) * t1576);
t1621 = t1624 * (pkin(2) * t1583 + pkin(3) * t1577);
t1620 = t1623 * (pkin(2) * t1584 + pkin(3) * t1578);
t1619 = t1625 * t1576;
t1618 = t1624 * t1577;
t1617 = t1623 * t1578;
t1567 = pkin(2) * t1585 + pkin(3) * t1579;
t1616 = t1567 * t1631;
t1615 = t1579 * t1628;
t1568 = pkin(2) * t1586 + pkin(3) * t1580;
t1614 = t1568 * t1630;
t1613 = t1580 * t1627;
t1569 = pkin(2) * t1587 + pkin(3) * t1581;
t1612 = t1569 * t1629;
t1611 = t1581 * t1626;
t1610 = t1567 * t1628;
t1609 = t1579 * t1631;
t1608 = t1568 * t1627;
t1607 = t1580 * t1630;
t1606 = t1569 * t1626;
t1605 = t1581 * t1629;
t1604 = 0.1e1 / pkin(2);
t1603 = 0.1e1 / pkin(3);
t1575 = t1596 * g(1) - t1593 * g(2);
t1574 = t1595 * g(1) - t1592 * g(2);
t1573 = t1594 * g(1) - t1591 * g(2);
t1572 = -t1593 * g(1) - t1596 * g(2);
t1571 = -t1592 * g(1) - t1595 * g(2);
t1570 = -t1591 * g(1) - t1594 * g(2);
t1554 = -g(3) * t1584 + t1575 * t1587;
t1553 = g(3) * t1587 + t1575 * t1584;
t1552 = -g(3) * t1583 + t1574 * t1586;
t1551 = g(3) * t1586 + t1574 * t1583;
t1550 = -g(3) * t1582 + t1573 * t1585;
t1549 = g(3) * t1585 + t1573 * t1582;
t1548 = -g(3) * t1578 + t1575 * t1581;
t1547 = g(3) * t1581 + t1575 * t1578;
t1546 = -g(3) * t1577 + t1574 * t1580;
t1545 = g(3) * t1580 + t1574 * t1577;
t1544 = -g(3) * t1576 + t1573 * t1579;
t1543 = g(3) * t1579 + t1573 * t1576;
t1 = [(t1591 * t1570 + t1592 * t1571 + t1593 * t1572) * MDP(1) - g(1) * MDP(8) + ((t1549 * t1615 + t1551 * t1613 + t1553 * t1611) * MDP(3) + (t1550 * t1615 + t1552 * t1613 + t1554 * t1611) * MDP(4) + (t1543 * t1615 + t1545 * t1613 + t1547 * t1611) * MDP(6) + (t1544 * t1615 + t1546 * t1613 + t1548 * t1611) * MDP(7) + ((-t1543 * t1610 - t1545 * t1608 - t1547 * t1606) * MDP(6) + (-t1544 * t1610 - t1546 * t1608 - t1548 * t1606) * MDP(7)) * t1603) * t1604; (t1594 * t1570 + t1595 * t1571 + t1596 * t1572) * MDP(1) - g(2) * MDP(8) + ((-t1549 * t1609 - t1551 * t1607 - t1553 * t1605) * MDP(3) + (-t1550 * t1609 - t1552 * t1607 - t1554 * t1605) * MDP(4) + (-t1543 * t1609 - t1545 * t1607 - t1547 * t1605) * MDP(6) + (-t1544 * t1609 - t1546 * t1607 - t1548 * t1605) * MDP(7) + ((t1543 * t1616 + t1545 * t1614 + t1547 * t1612) * MDP(6) + (t1544 * t1616 + t1546 * t1614 + t1548 * t1612) * MDP(7)) * t1603) * t1604; -g(3) * MDP(8) + ((-t1549 * t1619 - t1551 * t1618 - t1553 * t1617) * MDP(3) + (-t1550 * t1619 - t1552 * t1618 - t1554 * t1617) * MDP(4) + (-t1543 * t1619 - t1545 * t1618 - t1547 * t1617) * MDP(6) + (-t1544 * t1619 - t1546 * t1618 - t1548 * t1617) * MDP(7) + ((t1543 * t1622 + t1545 * t1621 + t1547 * t1620) * MDP(6) + (t1544 * t1622 + t1546 * t1621 + t1548 * t1620) * MDP(7)) * t1603) * t1604;];
taugX  = t1;
