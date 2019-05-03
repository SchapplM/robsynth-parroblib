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
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PRP1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:48
% EndTime: 2019-05-03 14:41:49
% DurationCPUTime: 0.52s
% Computational Cost: add. (654->143), mult. (1133->228), div. (36->3), fcn. (572->14), ass. (0->123)
t837 = 2 * pkin(2);
t798 = pkin(2) ^ 2;
t836 = 1 + t798;
t835 = mrSges(3,3) - mrSges(2,2);
t779 = legFrame(1,3);
t767 = sin(t779);
t834 = qJ(3,1) * t767;
t770 = cos(t779);
t833 = qJ(3,1) * t770;
t778 = legFrame(2,3);
t766 = sin(t778);
t832 = qJ(3,2) * t766;
t769 = cos(t778);
t831 = qJ(3,2) * t769;
t777 = legFrame(3,3);
t765 = sin(t777);
t830 = qJ(3,3) * t765;
t768 = cos(t777);
t829 = qJ(3,3) * t768;
t752 = -g(1) * t765 + g(2) * t768;
t755 = g(1) * t768 + g(2) * t765;
t758 = m(3) * qJ(3,3) + t835;
t761 = m(3) * pkin(2) + mrSges(2,1) + mrSges(3,1);
t780 = sin(qJ(2,3));
t783 = cos(qJ(2,3));
t716 = (-t752 * t761 - t755 * t758) * t783 - (t752 * t758 - t755 * t761) * t780;
t789 = qJ(3,3) ^ 2;
t762 = -t789 + t836;
t774 = t783 ^ 2;
t819 = t780 * t783;
t740 = 0.1e1 / (qJ(3,3) * t819 * t837 + t762 * t774 - t789 - t836);
t828 = t716 * t740;
t753 = -g(1) * t766 + g(2) * t769;
t756 = g(1) * t769 + g(2) * t766;
t759 = m(3) * qJ(3,2) + t835;
t781 = sin(qJ(2,2));
t784 = cos(qJ(2,2));
t717 = (-t753 * t761 - t756 * t759) * t784 - (t753 * t759 - t756 * t761) * t781;
t790 = qJ(3,2) ^ 2;
t763 = -t790 + t836;
t775 = t784 ^ 2;
t818 = t781 * t784;
t741 = 0.1e1 / (qJ(3,2) * t818 * t837 + t763 * t775 - t790 - t836);
t827 = t717 * t741;
t754 = -g(1) * t767 + g(2) * t770;
t757 = g(1) * t770 + g(2) * t767;
t760 = m(3) * qJ(3,1) + t835;
t782 = sin(qJ(2,1));
t785 = cos(qJ(2,1));
t718 = (-t754 * t761 - t757 * t760) * t785 - (t754 * t760 - t757 * t761) * t782;
t791 = qJ(3,1) ^ 2;
t764 = -t791 + t836;
t776 = t785 ^ 2;
t817 = t782 * t785;
t742 = 0.1e1 / (qJ(3,1) * t817 * t837 + t764 * t776 - t791 - t836);
t826 = t718 * t742;
t731 = -t752 * t783 + t755 * t780;
t825 = t731 * t740;
t732 = -t753 * t784 + t756 * t781;
t824 = t732 * t741;
t733 = -t754 * t785 + t757 * t782;
t823 = t733 * t742;
t822 = t740 * t752;
t821 = t741 * t753;
t820 = t742 * t754;
t816 = -0.2e1 * qJ(3,1) * t785;
t815 = -0.2e1 * qJ(3,2) * t784;
t814 = -0.2e1 * qJ(3,3) * t783;
t813 = pkin(2) * t834;
t812 = pkin(2) * t832;
t811 = pkin(2) * t830;
t810 = pkin(2) * t829;
t809 = pkin(2) * t831;
t808 = pkin(2) * t833;
t807 = -t764 * t767 + 0.2e1 * t808;
t806 = -t767 * t791 - t808;
t805 = t770 * t791 - t813;
t804 = -t763 * t766 + 0.2e1 * t809;
t803 = -t766 * t790 - t809;
t802 = t769 * t790 - t812;
t801 = -t762 * t765 + 0.2e1 * t810;
t800 = -t765 * t789 - t810;
t799 = t768 * t789 - t811;
t797 = koppelP(1,1);
t796 = koppelP(2,1);
t795 = koppelP(3,1);
t794 = koppelP(1,2);
t793 = koppelP(2,2);
t792 = koppelP(3,2);
t788 = mrSges(4,1);
t787 = mrSges(4,2);
t786 = xP(3);
t773 = m(1) + m(2) + m(3);
t772 = cos(t786);
t771 = sin(t786);
t751 = -t771 * t794 + t772 * t797;
t750 = -t771 * t793 + t772 * t796;
t749 = -t771 * t792 + t772 * t795;
t748 = -t771 * t797 - t772 * t794;
t747 = -t771 * t796 - t772 * t793;
t746 = -t771 * t795 - t772 * t792;
t745 = t764 * t770 + 0.2e1 * t813;
t744 = t763 * t769 + 0.2e1 * t812;
t743 = t762 * t768 + 0.2e1 * t811;
t739 = t770 * t816 + t782 * (pkin(2) * t770 + t834);
t738 = t767 * t816 + t782 * (pkin(2) * t767 - t833);
t737 = t769 * t815 + t781 * (pkin(2) * t769 + t832);
t736 = t766 * t815 + t781 * (pkin(2) * t766 - t831);
t735 = t768 * t814 + t780 * (pkin(2) * t768 + t830);
t734 = t765 * t814 + t780 * (pkin(2) * t765 - t829);
t730 = t805 * t785 - t782 * (t767 - t806);
t729 = t802 * t784 - t781 * (t766 - t803);
t728 = t799 * t783 - t780 * (t765 - t800);
t727 = t806 * t785 + t782 * (-t770 - t805);
t726 = t803 * t784 + t781 * (-t769 - t802);
t725 = t800 * t783 + t780 * (-t768 - t799);
t724 = -t745 * t817 + t798 * t767 + t807 * t776 + t767 - t808;
t723 = t745 * t776 - t770 * t798 + t807 * t817 - t770 - t813;
t722 = -t744 * t818 + t798 * t766 + t804 * t775 + t766 - t809;
t721 = t744 * t775 - t769 * t798 + t804 * t818 - t769 - t812;
t720 = -t743 * t819 + t798 * t765 + t801 * t774 + t765 - t810;
t719 = t743 * t774 - t768 * t798 + t801 * t819 - t768 - t811;
t1 = [t735 * t828 + t737 * t827 + t739 * t826 - g(1) * m(4) + (-t720 * t822 - t722 * t821 - t724 * t820) * t773 + (-t725 * t825 - t726 * t824 - t727 * t823) * m(3); t734 * t828 + t736 * t827 + t738 * t826 - g(2) * m(4) + (-t719 * t822 - t721 * t821 - t723 * t820) * t773 + (-t728 * t825 - t729 * t824 - t730 * t823) * m(3); -(-g(1) * t788 - g(2) * t787) * t771 + t772 * (g(1) * t787 - g(2) * t788) + (-(t723 * t751 + t724 * t748) * t754 * t773 + (t738 * t751 + t739 * t748) * t718 - (t727 * t748 + t730 * t751) * m(3) * t733) * t742 + (-(t721 * t750 + t722 * t747) * t753 * t773 + (t736 * t750 + t737 * t747) * t717 - (t726 * t747 + t729 * t750) * m(3) * t732) * t741 + (-(t719 * t749 + t720 * t746) * t752 * t773 + (t734 * t749 + t735 * t746) * t716 - (t725 * t746 + t728 * t749) * m(3) * t731) * t740;];
taugX  = t1;
