% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G2P2A0
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
%   see P3PRRR1G2P2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRR1G2P2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRRR1G2P2A0_gravload_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:37
% EndTime: 2020-03-09 21:18:38
% DurationCPUTime: 0.55s
% Computational Cost: add. (876->93), mult. (485->167), div. (126->5), fcn. (558->18), ass. (0->78)
t1615 = legFrame(3,2);
t1606 = sin(t1615);
t1612 = pkin(7) + qJ(2,3);
t1603 = qJ(3,3) + t1612;
t1588 = sin(t1603);
t1591 = cos(t1603);
t1597 = sin(t1612);
t1600 = cos(t1612);
t1640 = 0.1e1 / (-t1600 * t1588 + t1591 * t1597);
t1646 = t1606 * t1640;
t1616 = legFrame(2,2);
t1607 = sin(t1616);
t1613 = pkin(7) + qJ(2,2);
t1604 = qJ(3,2) + t1613;
t1589 = sin(t1604);
t1592 = cos(t1604);
t1598 = sin(t1613);
t1601 = cos(t1613);
t1639 = 0.1e1 / (-t1601 * t1589 + t1592 * t1598);
t1645 = t1607 * t1639;
t1617 = legFrame(1,2);
t1608 = sin(t1617);
t1614 = pkin(7) + qJ(2,1);
t1605 = qJ(3,1) + t1614;
t1590 = sin(t1605);
t1593 = cos(t1605);
t1599 = sin(t1614);
t1602 = cos(t1614);
t1638 = 0.1e1 / (-t1602 * t1590 + t1593 * t1599);
t1644 = t1608 * t1638;
t1609 = cos(t1615);
t1643 = t1609 * t1640;
t1610 = cos(t1616);
t1642 = t1610 * t1639;
t1611 = cos(t1617);
t1641 = t1611 * t1638;
t1637 = t1640 * (-pkin(2) * t1600 - pkin(3) * t1591);
t1636 = t1640 * t1591;
t1635 = t1639 * (-pkin(2) * t1601 - pkin(3) * t1592);
t1634 = t1639 * t1592;
t1633 = t1638 * (-pkin(2) * t1602 - pkin(3) * t1593);
t1632 = t1638 * t1593;
t1576 = pkin(2) * t1597 + pkin(3) * t1588;
t1631 = t1576 * t1643;
t1630 = t1588 * t1646;
t1577 = pkin(2) * t1598 + pkin(3) * t1589;
t1629 = t1577 * t1642;
t1628 = t1589 * t1645;
t1578 = pkin(2) * t1599 + pkin(3) * t1590;
t1627 = t1578 * t1641;
t1626 = t1590 * t1644;
t1625 = t1576 * t1646;
t1624 = t1588 * t1643;
t1623 = t1577 * t1645;
t1622 = t1589 * t1642;
t1621 = t1578 * t1644;
t1620 = t1590 * t1641;
t1619 = 0.1e1 / pkin(2);
t1618 = 0.1e1 / pkin(3);
t1587 = -t1608 * g(1) - t1611 * g(2);
t1586 = -t1611 * g(1) + t1608 * g(2);
t1585 = -t1607 * g(1) - t1610 * g(2);
t1584 = -t1610 * g(1) + t1607 * g(2);
t1583 = -t1606 * g(1) - t1609 * g(2);
t1582 = -t1609 * g(1) + t1606 * g(2);
t1566 = g(3) * t1602 - t1586 * t1599;
t1565 = g(3) * t1601 - t1584 * t1598;
t1564 = g(3) * t1600 - t1582 * t1597;
t1563 = g(3) * t1599 + t1586 * t1602;
t1562 = g(3) * t1598 + t1584 * t1601;
t1561 = g(3) * t1597 + t1582 * t1600;
t1560 = g(3) * t1593 - t1586 * t1590;
t1559 = g(3) * t1592 - t1584 * t1589;
t1558 = g(3) * t1591 - t1582 * t1588;
t1557 = g(3) * t1590 + t1586 * t1593;
t1556 = g(3) * t1589 + t1584 * t1592;
t1555 = g(3) * t1588 + t1582 * t1591;
t1 = [(t1606 * t1583 + t1607 * t1585 + t1608 * t1587) * MDP(1) - g(1) * MDP(8) + ((-t1561 * t1624 - t1562 * t1622 - t1563 * t1620) * MDP(3) + (-t1564 * t1624 - t1565 * t1622 - t1566 * t1620) * MDP(4) + (-t1555 * t1624 - t1556 * t1622 - t1557 * t1620) * MDP(6) + (-t1558 * t1624 - t1559 * t1622 - t1560 * t1620) * MDP(7) + ((t1555 * t1631 + t1556 * t1629 + t1557 * t1627) * MDP(6) + (t1558 * t1631 + t1559 * t1629 + t1560 * t1627) * MDP(7)) * t1618) * t1619; (t1609 * t1583 + t1610 * t1585 + t1611 * t1587) * MDP(1) - g(2) * MDP(8) + ((t1561 * t1630 + t1562 * t1628 + t1563 * t1626) * MDP(3) + (t1564 * t1630 + t1565 * t1628 + t1566 * t1626) * MDP(4) + (t1555 * t1630 + t1556 * t1628 + t1557 * t1626) * MDP(6) + (t1558 * t1630 + t1559 * t1628 + t1560 * t1626) * MDP(7) + ((-t1555 * t1625 - t1556 * t1623 - t1557 * t1621) * MDP(6) + (-t1558 * t1625 - t1559 * t1623 - t1560 * t1621) * MDP(7)) * t1618) * t1619; -g(3) * MDP(8) + ((-t1561 * t1636 - t1562 * t1634 - t1563 * t1632) * MDP(3) + (-t1564 * t1636 - t1565 * t1634 - t1566 * t1632) * MDP(4) + (-t1555 * t1636 - t1556 * t1634 - t1557 * t1632) * MDP(6) + (-t1558 * t1636 - t1559 * t1634 - t1560 * t1632) * MDP(7) + ((-t1555 * t1637 - t1556 * t1635 - t1557 * t1633) * MDP(6) + (-t1558 * t1637 - t1559 * t1635 - t1560 * t1633) * MDP(7)) * t1618) * t1619;];
taugX  = t1;
