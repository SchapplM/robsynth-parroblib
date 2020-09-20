% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRRR1G1P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [12x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRRRR1G1P3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRRRR1G1P3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1),zeros(12,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'P3PRRRR1G1P3A0_gravload_para_pf_mdp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:20
% EndTime: 2020-03-09 20:34:21
% DurationCPUTime: 0.50s
% Computational Cost: add. (169->73), mult. (381->154), div. (108->10), fcn. (450->18), ass. (0->82)
t1720 = legFrame(1,3);
t1702 = sin(t1720);
t1705 = cos(t1720);
t1699 = t1705 * g(1) + t1702 * g(2);
t1726 = sin(qJ(2,1));
t1732 = cos(qJ(2,1));
t1686 = -g(3) * t1732 + t1699 * t1726;
t1719 = legFrame(2,3);
t1701 = sin(t1719);
t1704 = cos(t1719);
t1698 = t1704 * g(1) + t1701 * g(2);
t1724 = sin(qJ(2,2));
t1730 = cos(qJ(2,2));
t1684 = -g(3) * t1730 + t1698 * t1724;
t1718 = legFrame(3,3);
t1700 = sin(t1718);
t1703 = cos(t1718);
t1697 = t1703 * g(1) + t1700 * g(2);
t1722 = sin(qJ(2,3));
t1728 = cos(qJ(2,3));
t1682 = -g(3) * t1728 + t1697 * t1722;
t1774 = MDP(1) * g(3);
t1733 = 0.1e1 / pkin(2);
t1770 = MDP(3) * t1733;
t1769 = MDP(4) * t1733;
t1768 = MDP(10) * t1733;
t1767 = MDP(11) * t1733;
t1694 = t1700 * g(1) - t1703 * g(2);
t1727 = cos(qJ(3,3));
t1712 = 0.1e1 / t1727;
t1721 = sin(qJ(3,3));
t1739 = g(3) * t1722 + t1697 * t1728;
t1766 = (t1721 * t1694 + t1727 * t1739) * t1712;
t1695 = t1701 * g(1) - t1704 * g(2);
t1729 = cos(qJ(3,2));
t1714 = 0.1e1 / t1729;
t1723 = sin(qJ(3,2));
t1738 = g(3) * t1724 + t1698 * t1730;
t1765 = (t1723 * t1695 + t1729 * t1738) * t1714;
t1696 = t1702 * g(1) - t1705 * g(2);
t1731 = cos(qJ(3,1));
t1716 = 0.1e1 / t1731;
t1725 = sin(qJ(3,1));
t1737 = g(3) * t1726 + t1699 * t1732;
t1764 = (t1725 * t1696 + t1731 * t1737) * t1716;
t1709 = 0.1e1 / t1722;
t1763 = t1682 * t1709;
t1710 = 0.1e1 / t1724;
t1762 = t1684 * t1710;
t1711 = 0.1e1 / t1726;
t1761 = t1686 * t1711;
t1757 = t1709 * t1712;
t1756 = t1709 / t1727 ^ 2;
t1755 = t1710 * t1714;
t1754 = t1710 / t1729 ^ 2;
t1753 = t1711 * t1716;
t1752 = t1711 / t1731 ^ 2;
t1751 = t1721 * t1728;
t1750 = t1723 * t1730;
t1749 = t1725 * t1732;
t1748 = t1727 * t1728;
t1747 = t1729 * t1730;
t1746 = t1731 * t1732;
t1688 = -t1700 * t1751 - t1703 * t1727;
t1745 = t1688 * t1756;
t1689 = -t1701 * t1750 - t1704 * t1729;
t1744 = t1689 * t1754;
t1690 = -t1702 * t1749 - t1705 * t1731;
t1743 = t1690 * t1752;
t1691 = -t1727 * t1700 + t1703 * t1751;
t1742 = t1691 * t1756;
t1692 = -t1729 * t1701 + t1704 * t1750;
t1741 = t1692 * t1754;
t1693 = -t1731 * t1702 + t1705 * t1749;
t1740 = t1693 * t1752;
t1736 = t1682 * t1721 * t1756;
t1735 = t1684 * t1723 * t1754;
t1734 = t1686 * t1725 * t1752;
t1674 = -t1731 * t1696 + t1725 * t1737;
t1672 = -t1729 * t1695 + t1723 * t1738;
t1670 = -t1727 * t1694 + t1721 * t1739;
t1 = [(-(t1702 * t1725 + t1705 * t1746) * t1753 - (t1701 * t1723 + t1704 * t1747) * t1755 - (t1700 * t1721 + t1703 * t1748) * t1757) * t1774 + (t1682 * t1745 + t1684 * t1744 + t1686 * t1743) * t1770 + (t1737 * t1743 + t1738 * t1744 + t1739 * t1745) * t1769 + ((t1674 * t1702 + t1690 * t1761) * t1716 + (t1672 * t1701 + t1689 * t1762) * t1714 + (t1670 * t1700 + t1688 * t1763) * t1712) * t1768 + (-t1688 * t1736 - t1689 * t1735 - t1690 * t1734 + t1700 * t1766 + t1701 * t1765 + t1702 * t1764) * t1767 - g(1) * MDP(12); (-(t1702 * t1746 - t1705 * t1725) * t1753 - (t1701 * t1747 - t1704 * t1723) * t1755 - (t1700 * t1748 - t1703 * t1721) * t1757) * t1774 + (t1682 * t1742 + t1684 * t1741 + t1686 * t1740) * t1770 + (t1737 * t1740 + t1738 * t1741 + t1739 * t1742) * t1769 + ((-t1674 * t1705 + t1693 * t1761) * t1716 + (-t1672 * t1704 + t1692 * t1762) * t1714 + (-t1670 * t1703 + t1691 * t1763) * t1712) * t1768 + (-t1691 * t1736 - t1692 * t1735 - t1693 * t1734 - t1703 * t1766 - t1704 * t1765 - t1705 * t1764) * t1767 - g(2) * MDP(12); (-0.3e1 * MDP(1) - MDP(12)) * g(3);];
taugX  = t1;
