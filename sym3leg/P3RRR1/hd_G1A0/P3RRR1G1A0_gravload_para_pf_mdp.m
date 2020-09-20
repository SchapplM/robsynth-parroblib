% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [10x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRR1G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(10,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'P3RRR1G1A0_gravload_para_pf_mdp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:47
% EndTime: 2019-05-03 15:38:48
% DurationCPUTime: 0.77s
% Computational Cost: add. (801->134), mult. (1004->249), div. (126->5), fcn. (1084->20), ass. (0->104)
t1777 = qJ(1,3) + qJ(2,3);
t1763 = sin(t1777);
t1766 = cos(t1777);
t1783 = sin(qJ(1,3));
t1786 = cos(qJ(1,3));
t1815 = 0.1e1 / (t1763 * t1786 - t1766 * t1783);
t1778 = qJ(1,2) + qJ(2,2);
t1764 = sin(t1778);
t1767 = cos(t1778);
t1784 = sin(qJ(1,2));
t1787 = cos(qJ(1,2));
t1814 = 0.1e1 / (t1764 * t1787 - t1767 * t1784);
t1779 = qJ(1,1) + qJ(2,1);
t1765 = sin(t1779);
t1768 = cos(t1779);
t1785 = sin(qJ(1,1));
t1788 = cos(qJ(1,1));
t1813 = 0.1e1 / (t1765 * t1788 - t1768 * t1785);
t1789 = xP(3);
t1775 = sin(t1789);
t1776 = cos(t1789);
t1790 = koppelP(3,2);
t1793 = koppelP(3,1);
t1740 = t1775 * t1793 + t1776 * t1790;
t1743 = -t1775 * t1790 + t1776 * t1793;
t1780 = legFrame(3,3);
t1769 = sin(t1780);
t1772 = cos(t1780);
t1698 = (t1740 * t1772 - t1743 * t1769) * t1766 - (t1740 * t1769 + t1743 * t1772) * t1763;
t1812 = t1698 * t1815;
t1791 = koppelP(2,2);
t1794 = koppelP(2,1);
t1741 = t1775 * t1794 + t1776 * t1791;
t1744 = -t1775 * t1791 + t1776 * t1794;
t1781 = legFrame(2,3);
t1770 = sin(t1781);
t1773 = cos(t1781);
t1699 = (t1741 * t1773 - t1744 * t1770) * t1767 - (t1741 * t1770 + t1744 * t1773) * t1764;
t1811 = t1699 * t1814;
t1792 = koppelP(1,2);
t1795 = koppelP(1,1);
t1742 = t1775 * t1795 + t1776 * t1792;
t1745 = -t1775 * t1792 + t1776 * t1795;
t1782 = legFrame(1,3);
t1771 = sin(t1782);
t1774 = cos(t1782);
t1700 = (t1742 * t1774 - t1745 * t1771) * t1768 - (t1742 * t1771 + t1745 * t1774) * t1765;
t1810 = t1700 * t1813;
t1746 = -t1769 * g(1) + g(2) * t1772;
t1749 = g(1) * t1772 + g(2) * t1769;
t1713 = -t1746 * t1766 + t1749 * t1763;
t1809 = t1713 * t1815;
t1714 = t1746 * t1763 + t1749 * t1766;
t1808 = t1714 * t1815;
t1747 = -t1770 * g(1) + g(2) * t1773;
t1750 = g(1) * t1773 + g(2) * t1770;
t1715 = -t1747 * t1767 + t1750 * t1764;
t1807 = t1715 * t1814;
t1716 = t1747 * t1764 + t1750 * t1767;
t1806 = t1716 * t1814;
t1748 = -t1771 * g(1) + g(2) * t1774;
t1751 = g(1) * t1774 + g(2) * t1771;
t1717 = -t1748 * t1768 + t1751 * t1765;
t1805 = t1717 * t1813;
t1718 = t1748 * t1765 + t1751 * t1768;
t1804 = t1718 * t1813;
t1725 = t1763 * t1772 + t1766 * t1769;
t1803 = t1725 * t1815;
t1726 = t1763 * t1769 - t1766 * t1772;
t1802 = t1726 * t1815;
t1727 = t1764 * t1773 + t1767 * t1770;
t1801 = t1727 * t1814;
t1728 = t1764 * t1770 - t1767 * t1773;
t1800 = t1728 * t1814;
t1729 = t1765 * t1774 + t1768 * t1771;
t1799 = t1729 * t1813;
t1730 = t1765 * t1771 - t1768 * t1774;
t1798 = t1730 * t1813;
t1797 = 0.1e1 / pkin(1);
t1796 = 0.1e1 / pkin(2);
t1759 = t1785 * t1792 + t1788 * t1795;
t1758 = t1785 * t1795 - t1788 * t1792;
t1757 = t1784 * t1791 + t1787 * t1794;
t1756 = t1784 * t1794 - t1787 * t1791;
t1755 = t1783 * t1790 + t1786 * t1793;
t1754 = t1783 * t1793 - t1786 * t1790;
t1753 = g(1) * t1776 + g(2) * t1775;
t1752 = g(1) * t1775 - g(2) * t1776;
t1724 = t1748 * t1785 + t1751 * t1788;
t1723 = -t1748 * t1788 + t1751 * t1785;
t1722 = t1747 * t1784 + t1750 * t1787;
t1721 = -t1747 * t1787 + t1750 * t1784;
t1720 = t1746 * t1783 + t1749 * t1786;
t1719 = -t1746 * t1786 + t1749 * t1783;
t1706 = pkin(1) * (-t1771 * t1785 + t1774 * t1788) - t1730 * pkin(2);
t1705 = pkin(1) * (-t1770 * t1784 + t1773 * t1787) - t1728 * pkin(2);
t1704 = pkin(1) * (-t1769 * t1783 + t1772 * t1786) - t1726 * pkin(2);
t1703 = pkin(1) * (t1771 * t1788 + t1774 * t1785) + t1729 * pkin(2);
t1702 = pkin(1) * (t1770 * t1787 + t1773 * t1784) + t1727 * pkin(2);
t1701 = pkin(1) * (t1769 * t1786 + t1772 * t1783) + t1725 * pkin(2);
t1697 = pkin(1) * ((t1758 * t1776 - t1759 * t1775) * t1774 + t1771 * (t1758 * t1775 + t1759 * t1776)) - t1700 * pkin(2);
t1696 = pkin(1) * ((t1756 * t1776 - t1757 * t1775) * t1773 + t1770 * (t1756 * t1775 + t1757 * t1776)) - t1699 * pkin(2);
t1695 = pkin(1) * ((t1754 * t1776 - t1755 * t1775) * t1772 + t1769 * (t1754 * t1775 + t1755 * t1776)) - t1698 * pkin(2);
t1 = [(-t1752 * t1775 - t1753 * t1776) * MDP(10) + ((-t1719 * t1802 - t1721 * t1800 - t1723 * t1798) * MDP(2) + (-t1720 * t1802 - t1722 * t1800 - t1724 * t1798) * MDP(3) + (-t1713 * t1802 - t1715 * t1800 - t1717 * t1798) * MDP(5) + (-t1714 * t1802 - t1716 * t1800 - t1718 * t1798) * MDP(6) + ((-t1704 * t1809 - t1705 * t1807 - t1706 * t1805) * MDP(5) + (-t1704 * t1808 - t1705 * t1806 - t1706 * t1804) * MDP(6)) * t1796) * t1797; (t1752 * t1776 - t1753 * t1775) * MDP(10) + ((t1719 * t1803 + t1721 * t1801 + t1723 * t1799) * MDP(2) + (t1720 * t1803 + t1722 * t1801 + t1724 * t1799) * MDP(3) + (t1713 * t1803 + t1715 * t1801 + t1717 * t1799) * MDP(5) + (t1714 * t1803 + t1716 * t1801 + t1718 * t1799) * MDP(6) + ((-t1701 * t1809 - t1702 * t1807 - t1703 * t1805) * MDP(5) + (-t1701 * t1808 - t1702 * t1806 - t1703 * t1804) * MDP(6)) * t1796) * t1797; t1752 * MDP(8) + t1753 * MDP(9) + ((-t1719 * t1812 - t1721 * t1811 - t1723 * t1810) * MDP(2) + (-t1720 * t1812 - t1722 * t1811 - t1724 * t1810) * MDP(3) + (-t1698 * t1809 - t1699 * t1807 - t1700 * t1805) * MDP(5) + (-t1698 * t1808 - t1699 * t1806 - t1700 * t1804) * MDP(6) + ((-t1695 * t1809 - t1696 * t1807 - t1697 * t1805) * MDP(5) + (-t1695 * t1808 - t1696 * t1806 - t1697 * t1804) * MDP(6)) * t1796) * t1797;];
taugX  = t1;
