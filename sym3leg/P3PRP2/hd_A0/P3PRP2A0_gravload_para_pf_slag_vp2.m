% Calculate Gravitation load for parallel robot
% P3PRP2A0
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
%   mass of all robot links (including platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PRP2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:39:05
% EndTime: 2018-12-20 17:39:05
% DurationCPUTime: 0.48s
% Computational Cost: add. (654->143), mult. (1133->233), div. (36->3), fcn. (572->14), ass. (0->125)
t811 = (pkin(2) ^ 2);
t846 = -t811 - 1;
t845 = mrSges(3,3) - mrSges(2,2);
t792 = legFrame(1,3);
t780 = sin(t792);
t844 = qJ(3,1) * t780;
t783 = cos(t792);
t843 = qJ(3,1) * t783;
t791 = legFrame(2,3);
t779 = sin(t791);
t842 = qJ(3,2) * t779;
t782 = cos(t791);
t841 = qJ(3,2) * t782;
t790 = legFrame(3,3);
t778 = sin(t790);
t840 = qJ(3,3) * t778;
t781 = cos(t790);
t839 = qJ(3,3) * t781;
t759 = -t778 * g(1) + t781 * g(2);
t762 = t781 * g(1) + t778 * g(2);
t771 = m(3) * qJ(3,3) + t845;
t774 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t793 = sin(qJ(2,3));
t796 = cos(qJ(2,3));
t720 = (-t759 * t771 - t774 * t762) * t796 + (t759 * t774 - t771 * t762) * t793;
t802 = qJ(3,3) ^ 2;
t775 = -t802 - t846;
t787 = t796 ^ 2;
t826 = t796 * qJ(3,3);
t820 = 0.2e1 * t826;
t744 = 0.1e1 / (pkin(2) * t793 * t820 + t775 * t787 - t802 + t846);
t838 = t720 * t744;
t760 = -t779 * g(1) + t782 * g(2);
t763 = t782 * g(1) + t779 * g(2);
t772 = m(3) * qJ(3,2) + t845;
t794 = sin(qJ(2,2));
t797 = cos(qJ(2,2));
t721 = (-t760 * t772 - t774 * t763) * t797 + (t760 * t774 - t772 * t763) * t794;
t803 = qJ(3,2) ^ 2;
t776 = -t803 - t846;
t788 = t797 ^ 2;
t825 = t797 * qJ(3,2);
t819 = 0.2e1 * t825;
t745 = 0.1e1 / (pkin(2) * t794 * t819 + t776 * t788 - t803 + t846);
t837 = t721 * t745;
t761 = -t780 * g(1) + t783 * g(2);
t764 = t783 * g(1) + t780 * g(2);
t773 = m(3) * qJ(3,1) + t845;
t795 = sin(qJ(2,1));
t798 = cos(qJ(2,1));
t722 = (-t761 * t773 - t774 * t764) * t798 + (t761 * t774 - t773 * t764) * t795;
t804 = qJ(3,1) ^ 2;
t777 = -t804 - t846;
t789 = t798 ^ 2;
t824 = t798 * qJ(3,1);
t818 = 0.2e1 * t824;
t746 = 0.1e1 / (pkin(2) * t795 * t818 + t777 * t789 - t804 + t846);
t836 = t722 * t746;
t735 = -t793 * t759 + t796 * t762;
t835 = t735 * t744;
t736 = -t794 * t760 + t797 * t763;
t834 = t736 * t745;
t737 = -t795 * t761 + t798 * t764;
t833 = t737 * t746;
t832 = t744 * t762;
t831 = t745 * t763;
t830 = t746 * t764;
t829 = t793 * t796;
t828 = t794 * t797;
t827 = t795 * t798;
t765 = pkin(2) * t840;
t823 = t802 * t781 + t765;
t766 = pkin(2) * t842;
t822 = t803 * t782 + t766;
t767 = pkin(2) * t844;
t821 = t804 * t783 + t767;
t817 = pkin(2) * t843;
t816 = pkin(2) * t841;
t815 = pkin(2) * t839;
t814 = t780 * t804 - t817;
t813 = t779 * t803 - t816;
t812 = t778 * t802 - t815;
t810 = koppelP(1,1);
t809 = koppelP(2,1);
t808 = koppelP(3,1);
t807 = koppelP(1,2);
t806 = koppelP(2,2);
t805 = koppelP(3,2);
t801 = mrSges(4,1);
t800 = mrSges(4,2);
t799 = xP(3);
t786 = m(1) + m(2) + m(3);
t785 = cos(t799);
t784 = sin(t799);
t758 = -t784 * t807 + t785 * t810;
t757 = -t784 * t806 + t785 * t809;
t756 = -t784 * t805 + t785 * t808;
t755 = -t784 * t810 - t785 * t807;
t754 = -t784 * t809 - t785 * t806;
t753 = -t784 * t808 - t785 * t805;
t752 = t777 * t780 + 0.2e1 * t817;
t751 = t776 * t779 + 0.2e1 * t816;
t750 = t775 * t778 + 0.2e1 * t815;
t749 = t783 * t777 - 0.2e1 * t767;
t748 = t782 * t776 - 0.2e1 * t766;
t747 = t781 * t775 - 0.2e1 * t765;
t743 = -0.2e1 * t783 * t824 + t795 * (pkin(2) * t783 - t844);
t742 = t780 * t818 - t795 * (pkin(2) * t780 + t843);
t741 = -0.2e1 * t782 * t825 + t794 * (pkin(2) * t782 - t842);
t740 = t779 * t819 - t794 * (pkin(2) * t779 + t841);
t739 = -0.2e1 * t781 * t826 + t793 * (pkin(2) * t781 - t840);
t738 = t778 * t820 - t793 * (pkin(2) * t778 + t839);
t734 = t814 * t798 - t795 * (t783 + t821);
t733 = t813 * t797 - t794 * (t782 + t822);
t732 = t812 * t796 - t793 * (t781 + t823);
t731 = t821 * t798 - t795 * (-t780 - t814);
t730 = t822 * t797 - t794 * (-t779 - t813);
t729 = t823 * t796 - t793 * (-t778 - t812);
t728 = -t749 * t827 + t752 * t789 - t780 * t811 - t780 - t817;
t727 = -t748 * t828 + t751 * t788 - t779 * t811 - t779 - t816;
t726 = -t747 * t829 + t750 * t787 - t778 * t811 - t778 - t815;
t725 = t749 * t789 + t752 * t827 - t811 * t783 + t767 - t783;
t724 = t748 * t788 + t751 * t828 - t811 * t782 + t766 - t782;
t723 = t747 * t787 + t750 * t829 - t811 * t781 + t765 - t781;
t1 = [t738 * t838 + t740 * t837 + t742 * t836 - g(1) * m(4) + (-t723 * t832 - t724 * t831 - t725 * t830) * t786 + (t729 * t835 + t730 * t834 + t731 * t833) * m(3); t739 * t838 + t741 * t837 + t743 * t836 - g(2) * m(4) + (-t726 * t832 - t727 * t831 - t728 * t830) * t786 + (t732 * t835 + t733 * t834 + t734 * t833) * m(3); -(-g(1) * t801 - g(2) * t800) * t784 + t785 * (g(1) * t800 - g(2) * t801) + (-(t725 * t755 + t728 * t758) * t764 * t786 + (t742 * t755 + t743 * t758) * t722 + (t731 * t755 + t734 * t758) * m(3) * t737) * t746 + (-(t724 * t754 + t727 * t757) * t763 * t786 + (t740 * t754 + t741 * t757) * t721 + (t730 * t754 + t733 * t757) * m(3) * t736) * t745 + (-(t723 * t753 + t726 * t756) * t762 * t786 + (t738 * t753 + t739 * t756) * t720 + (t729 * t753 + t732 * t756) * m(3) * t735) * t744;];
taugX  = t1;
