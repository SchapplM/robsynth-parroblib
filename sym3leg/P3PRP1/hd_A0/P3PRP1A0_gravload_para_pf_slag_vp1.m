% Calculate Gravitation load for parallel robot
% P3PRP1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRP1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:38
% EndTime: 2019-05-03 14:41:39
% DurationCPUTime: 0.56s
% Computational Cost: add. (690->153), mult. (1242->252), div. (36->3), fcn. (644->14), ass. (0->125)
t852 = 2 * pkin(2);
t811 = pkin(2) ^ 2;
t851 = 1 + t811;
t791 = legFrame(1,3);
t776 = sin(t791);
t850 = qJ(3,1) * t776;
t779 = cos(t791);
t849 = qJ(3,1) * t779;
t790 = legFrame(2,3);
t775 = sin(t790);
t848 = qJ(3,2) * t775;
t778 = cos(t790);
t847 = qJ(3,2) * t778;
t789 = legFrame(3,3);
t774 = sin(t789);
t846 = qJ(3,3) * t774;
t777 = cos(t789);
t845 = qJ(3,3) * t777;
t765 = -t774 * g(1) + t777 * g(2);
t768 = t777 * g(1) + t774 * g(2);
t786 = rSges(3,3) + qJ(3,3);
t792 = sin(qJ(2,3));
t795 = cos(qJ(2,3));
t798 = pkin(2) + rSges(3,1);
t729 = ((-t765 * t798 - t786 * t768) * m(3) + m(2) * (-rSges(2,1) * t765 + rSges(2,2) * t768)) * t795 + t792 * ((-t765 * t786 + t798 * t768) * m(3) + m(2) * (rSges(2,1) * t768 + rSges(2,2) * t765));
t802 = qJ(3,3) ^ 2;
t771 = -t802 + t851;
t783 = t795 ^ 2;
t832 = t795 * qJ(3,3);
t753 = 0.1e1 / (t792 * t832 * t852 + t771 * t783 - t802 - t851);
t844 = t729 * t753;
t766 = -t775 * g(1) + t778 * g(2);
t769 = t778 * g(1) + t775 * g(2);
t787 = rSges(3,3) + qJ(3,2);
t793 = sin(qJ(2,2));
t796 = cos(qJ(2,2));
t730 = ((-t766 * t798 - t787 * t769) * m(3) + m(2) * (-rSges(2,1) * t766 + rSges(2,2) * t769)) * t796 + t793 * ((-t766 * t787 + t798 * t769) * m(3) + m(2) * (rSges(2,1) * t769 + rSges(2,2) * t766));
t803 = qJ(3,2) ^ 2;
t772 = -t803 + t851;
t784 = t796 ^ 2;
t831 = t796 * qJ(3,2);
t754 = 0.1e1 / (t793 * t831 * t852 + t772 * t784 - t803 - t851);
t843 = t730 * t754;
t767 = -t776 * g(1) + t779 * g(2);
t770 = t779 * g(1) + t776 * g(2);
t788 = rSges(3,3) + qJ(3,1);
t794 = sin(qJ(2,1));
t797 = cos(qJ(2,1));
t731 = ((-t767 * t798 - t788 * t770) * m(3) + m(2) * (-rSges(2,1) * t767 + rSges(2,2) * t770)) * t797 + t794 * ((-t767 * t788 + t798 * t770) * m(3) + m(2) * (rSges(2,1) * t770 + rSges(2,2) * t767));
t804 = qJ(3,1) ^ 2;
t773 = -t804 + t851;
t785 = t797 ^ 2;
t830 = t797 * qJ(3,1);
t755 = 0.1e1 / (t794 * t830 * t852 + t773 * t785 - t804 - t851);
t842 = t731 * t755;
t744 = -t765 * t795 + t768 * t792;
t841 = t744 * t753;
t745 = -t766 * t796 + t769 * t793;
t840 = t745 * t754;
t746 = -t767 * t797 + t770 * t794;
t839 = t746 * t755;
t838 = t753 * t765;
t837 = t754 * t766;
t836 = t755 * t767;
t835 = t792 * t795;
t834 = t793 * t796;
t833 = t794 * t797;
t829 = -0.2e1 * t832;
t828 = -0.2e1 * t831;
t827 = -0.2e1 * t830;
t826 = pkin(2) * t850;
t825 = pkin(2) * t848;
t824 = pkin(2) * t846;
t823 = pkin(2) * t845;
t822 = pkin(2) * t847;
t821 = pkin(2) * t849;
t820 = -t776 * t773 + 0.2e1 * t821;
t819 = -t804 * t776 - t821;
t818 = t779 * t804 - t826;
t817 = -t775 * t772 + 0.2e1 * t822;
t816 = -t803 * t775 - t822;
t815 = t778 * t803 - t825;
t814 = -t774 * t771 + 0.2e1 * t823;
t813 = -t802 * t774 - t823;
t812 = t777 * t802 - t824;
t810 = koppelP(1,1);
t809 = koppelP(2,1);
t808 = koppelP(3,1);
t807 = koppelP(1,2);
t806 = koppelP(2,2);
t805 = koppelP(3,2);
t801 = rSges(4,1);
t800 = rSges(4,2);
t799 = xP(3);
t782 = m(1) + m(2) + m(3);
t781 = cos(t799);
t780 = sin(t799);
t764 = -t780 * t807 + t781 * t810;
t763 = -t780 * t806 + t781 * t809;
t762 = -t780 * t805 + t781 * t808;
t761 = -t780 * t810 - t781 * t807;
t760 = -t780 * t809 - t781 * t806;
t759 = -t780 * t808 - t781 * t805;
t758 = t773 * t779 + 0.2e1 * t826;
t757 = t772 * t778 + 0.2e1 * t825;
t756 = t771 * t777 + 0.2e1 * t824;
t752 = t779 * t827 + t794 * (pkin(2) * t779 + t850);
t751 = t776 * t827 + t794 * (pkin(2) * t776 - t849);
t750 = t778 * t828 + t793 * (pkin(2) * t778 + t848);
t749 = t775 * t828 + t793 * (pkin(2) * t775 - t847);
t748 = t777 * t829 + t792 * (pkin(2) * t777 + t846);
t747 = t774 * t829 + t792 * (pkin(2) * t774 - t845);
t743 = t818 * t797 - t794 * (t776 - t819);
t742 = t815 * t796 - t793 * (t775 - t816);
t741 = t812 * t795 - t792 * (t774 - t813);
t740 = t819 * t797 + t794 * (-t779 - t818);
t739 = t816 * t796 + t793 * (-t778 - t815);
t738 = t813 * t795 + t792 * (-t777 - t812);
t737 = -t758 * t833 + t811 * t776 + t820 * t785 + t776 - t821;
t736 = t758 * t785 - t779 * t811 + t820 * t833 - t779 - t826;
t735 = -t757 * t834 + t811 * t775 + t817 * t784 + t775 - t822;
t734 = t757 * t784 - t778 * t811 + t817 * t834 - t778 - t825;
t733 = -t756 * t835 + t811 * t774 + t814 * t783 + t774 - t823;
t732 = t756 * t783 - t777 * t811 + t814 * t835 - t777 - t824;
t1 = [t748 * t844 + t750 * t843 + t752 * t842 - m(4) * g(1) + (-t733 * t838 - t735 * t837 - t737 * t836) * t782 + (-t738 * t841 - t739 * t840 - t740 * t839) * m(3); t747 * t844 + t749 * t843 + t751 * t842 - m(4) * g(2) + (-t732 * t838 - t734 * t837 - t736 * t836) * t782 + (-t741 * t841 - t742 * t840 - t743 * t839) * m(3); m(4) * ((g(1) * t801 + g(2) * t800) * t780 + (g(1) * t800 - g(2) * t801) * t781) + (-(t736 * t764 + t737 * t761) * t767 * t782 + (t751 * t764 + t752 * t761) * t731 - (t740 * t761 + t743 * t764) * m(3) * t746) * t755 + (-(t734 * t763 + t735 * t760) * t766 * t782 + (t749 * t763 + t750 * t760) * t730 - (t739 * t760 + t742 * t763) * m(3) * t745) * t754 + (-(t732 * t762 + t733 * t759) * t765 * t782 + (t747 * t762 + t748 * t759) * t729 - (t738 * t759 + t741 * t762) * m(3) * t744) * t753;];
taugX  = t1;
