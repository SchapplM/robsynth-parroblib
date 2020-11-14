% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR9V1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:07
% EndTime: 2020-08-06 18:47:12
% DurationCPUTime: 5.68s
% Computational Cost: add. (17001->432), mult. (14874->677), div. (3786->11), fcn. (12975->53), ass. (0->295)
t753 = 2 * pkin(7);
t705 = t753 + qJ(3,1);
t737 = cos(qJ(3,1));
t928 = (cos(t705) + t737) * pkin(2);
t704 = t753 + qJ(3,2);
t735 = cos(qJ(3,2));
t927 = (cos(t704) + t735) * pkin(2);
t703 = t753 + qJ(3,3);
t733 = cos(qJ(3,3));
t926 = (cos(t703) + t733) * pkin(2);
t708 = pkin(7) + qJ(3,1);
t659 = 2 * t708;
t925 = (cos(t659) + 0.1e1) * pkin(3);
t707 = pkin(7) + qJ(3,2);
t658 = 2 * t707;
t924 = (cos(t658) + 0.1e1) * pkin(3);
t706 = pkin(7) + qJ(3,3);
t657 = 2 * t706;
t923 = (cos(t657) + 0.1e1) * pkin(3);
t722 = pkin(5) + qJ(2,3);
t696 = -pkin(6) - t722;
t688 = 0.1e1 / t696;
t739 = xDP(3);
t644 = sin(t657);
t666 = sin(t706);
t727 = sin(qJ(3,3));
t906 = 2 * pkin(1);
t623 = t666 * t906 + pkin(3) * t644 + (sin(t703) + t727) * pkin(2);
t672 = cos(t706);
t660 = 0.1e1 / t672;
t802 = t623 * t660 / 0.2e1;
t718 = cos(pkin(7));
t675 = t718 * pkin(2);
t656 = t675 + pkin(1);
t728 = sin(qJ(1,3));
t734 = cos(qJ(1,3));
t626 = t728 * t656 + t734 * t696;
t629 = t656 * t734 - t728 * t696;
t719 = legFrame(3,3);
t676 = sin(t719);
t679 = cos(t719);
t635 = t676 * t734 + t679 * t728;
t896 = pkin(3) * t672;
t608 = t626 * t679 + t629 * t676 + t635 * t896;
t740 = xDP(2);
t874 = t608 * t740;
t632 = -t676 * t728 + t679 * t734;
t605 = -t626 * t676 + t629 * t679 + t632 * t896;
t741 = xDP(1);
t877 = t605 * t741;
t596 = (t739 * t802 + t874 + t877) * t688;
t864 = t660 * t739;
t614 = (t632 * t741 + t635 * t740 + t666 * t864) * t688;
t612 = pkin(1) * t614;
t752 = 3 * pkin(7);
t756 = 2 * qJ(3,3);
t759 = pkin(3) ^ 2;
t762 = pkin(2) ^ 2;
t760 = 0.1e1 / pkin(3);
t845 = t739 * t760;
t790 = pkin(3) * t845 / 0.2e1;
t763 = pkin(1) ^ 2;
t803 = -0.3e1 * t759 - 0.2e1 * t762 - (4 * t763);
t804 = -0.4e1 * pkin(3) * t675;
t814 = t623 * t864;
t843 = -0.2e1 * pkin(2) * pkin(3);
t865 = t660 * t688;
t871 = t614 * t688;
t897 = -t739 / 0.2e1;
t905 = 0.1e1 / t672 ^ 2;
t716 = t739 ^ 2;
t910 = t716 * t760;
t913 = t644 * t696;
t920 = t926 / 0.2e1 + t923 / 0.2e1;
t581 = (t656 + t896) * t905 * t910 * t865 + (t905 * t897 * t913 + (-(-t759 * cos((3 * t706)) + (-0.4e1 * t696 ^ 2 + t803) * t672 + t804 + (cos((t756 + t752)) + cos((t756 + pkin(7)))) * t843 + (-cos((-pkin(7) + qJ(3,3))) - cos((t752 + qJ(3,3)))) * t762) * t614 / 0.4e1 - t660 * t790 * t913 + (-t923 - t926) * (-t612 - (-t877 / 0.2e1 - t874 / 0.2e1 - t814 / 0.4e1) * t688) - (t672 * t906 + t920) * t596) * t660) * t871;
t661 = t660 * t905;
t772 = t814 + 0.2e1 * t874 + 0.2e1 * t877;
t799 = -t865 / 0.2e1;
t587 = (-t661 * t910 - (-t596 + t660 * t614 * t920 + (t612 * t660 + t772 * t799) * t672) * t614) * t688;
t725 = Ifges(3,1) - Ifges(3,2);
t851 = t725 * t727;
t602 = t614 * t851;
t650 = t722 * mrSges(3,2) - Ifges(3,6);
t653 = t722 * mrSges(3,1) - Ifges(3,5);
t717 = sin(pkin(7));
t617 = -(-t650 * t733 - t653 * t727) * t718 - t717 * (t650 * t727 - t653 * t733);
t742 = (m(2) + m(3));
t695 = pkin(1) * t742;
t890 = mrSges(3,2) * t727;
t779 = mrSges(3,1) * t733 - t890;
t682 = t727 * mrSges(3,1);
t685 = mrSges(3,2) * t733;
t839 = t685 + t682;
t881 = m(3) * pkin(2) + mrSges(2,1);
t620 = -(-t779 - t881) * t718 - (mrSges(2,2) + t839) * t717 + t695;
t726 = mrSges(2,3) + mrSges(3,3);
t818 = m(3) * pkin(5) + t726;
t638 = (t742 * qJ(2,3)) + t818;
t898 = pkin(2) * mrSges(3,2);
t693 = -t898 / 0.2e1;
t899 = pkin(2) * mrSges(3,1);
t694 = t899 / 0.4e1;
t702 = t718 ^ 2;
t710 = 0.2e1 * t899;
t713 = t733 ^ 2;
t743 = Ifges(3,6) / 0.2e1;
t744 = Ifges(3,5) / 0.4e1;
t754 = -pkin(5) / 0.4e1;
t755 = -pkin(5) / 0.2e1;
t776 = -0.2e1 * mrSges(3,3) * pkin(5) - (t742 * t763) - Ifges(2,1) - Ifges(3,2) - Ifges(1,3);
t780 = t762 * m(3) - Ifges(2,1) + Ifges(2,2) + t725;
t854 = t725 * t713;
t870 = t614 * t717;
t784 = t854 * t870;
t855 = t716 / pkin(3) ^ 2;
t787 = t661 * t666 * t855;
t901 = -0.2e1 * t717;
t791 = mrSges(3,1) * t906 + t898 * t901;
t795 = 0.4e1 * t718 * t845;
t805 = t717 * t845;
t810 = t660 * t845;
t824 = pkin(2) * t890;
t817 = (t824 - t725) * t870;
t820 = -t899 / 0.2e1;
t900 = mrSges(3,2) * pkin(1);
t821 = -t900 / 0.2e1;
t827 = mrSges(3,1) * t612;
t887 = Ifges(3,4) * t727;
t833 = 0.4e1 * t887;
t834 = 0.4e1 * Ifges(3,4) * t717;
t835 = t717 * t906;
t836 = -0.2e1 * t900;
t848 = t727 * t733;
t884 = t713 * Ifges(3,4);
t893 = Ifges(2,4) - Ifges(3,4);
t914 = 0.2e1 * t884 - Ifges(3,4);
t842 = -(-((t710 + t833) * t733 - 0.2e1 * t824 + t780) * t702 - (t713 * t834 + t791 * t733 + (t881 - t890) * t906) * t718 - t854 + (t682 + mrSges(2,2)) * t835 - t722 ^ 2 * m(3) - (m(2) * qJ(2,3) ^ 2) + t776) * t587 + t620 * t581 - t617 * t787 - 0.2e1 * (Ifges(3,4) * t848 + t702 * t854 - t726 * qJ(2,3) + (-(-pkin(2) * t682 + t725 * t848 + t893) * t718 + pkin(1) * t685) * t717) * t587 - 0.2e1 * ((((-qJ(2,3) / 0.2e1 + t755) * mrSges(3,2) + t743) * t810 - t827) * t717 - t602) * t733 * t810 - 0.2e1 * t614 * (Ifges(3,4) * t810 - t638 * t596) - 0.4e1 * (((t693 + t851) * t733 + t727 * t820 + t914) * t702 - t884) * t614 * t810 + ((-t784 + (((-qJ(2,3) / 0.4e1 + t754) * mrSges(3,1) + t744) * t810 - ((t694 + t887) * t901 + t821) * t614) * t733 - t817 / 0.2e1 - t727 * (-t650 * t810 - 0.2e1 * t827) / 0.4e1) * t795 - (-t614 * t836 - t653 * t810) * t727 * t805) * t660;
t723 = pkin(5) + qJ(2,2);
t697 = -pkin(6) - t723;
t689 = 0.1e1 / t697;
t645 = sin(t658);
t667 = sin(t707);
t729 = sin(qJ(3,2));
t624 = t667 * t906 + pkin(3) * t645 + (sin(t704) + t729) * pkin(2);
t673 = cos(t707);
t662 = 0.1e1 / t673;
t801 = t624 * t662 / 0.2e1;
t730 = sin(qJ(1,2));
t736 = cos(qJ(1,2));
t627 = t730 * t656 + t736 * t697;
t630 = t656 * t736 - t730 * t697;
t720 = legFrame(2,3);
t677 = sin(t720);
t680 = cos(t720);
t636 = t677 * t736 + t680 * t730;
t895 = pkin(3) * t673;
t609 = t627 * t680 + t630 * t677 + t636 * t895;
t873 = t609 * t740;
t633 = -t677 * t730 + t680 * t736;
t606 = -t627 * t677 + t630 * t680 + t633 * t895;
t876 = t606 * t741;
t597 = (t739 * t801 + t873 + t876) * t689;
t861 = t662 * t739;
t615 = (t633 * t741 + t636 * t740 + t667 * t861) * t689;
t613 = pkin(1) * t615;
t757 = 2 * qJ(3,2);
t813 = t624 * t861;
t862 = t662 * t689;
t869 = t615 * t689;
t904 = 0.1e1 / t673 ^ 2;
t912 = t645 * t697;
t921 = t927 / 0.2e1 + t924 / 0.2e1;
t582 = (t656 + t895) * t904 * t910 * t862 + (t904 * t897 * t912 + (-(-t759 * cos((3 * t707)) + (-0.4e1 * t697 ^ 2 + t803) * t673 + t804 + (cos((t752 + t757)) + cos((t757 + pkin(7)))) * t843 + (-cos((qJ(3,2) + t752)) - cos((-pkin(7) + qJ(3,2)))) * t762) * t615 / 0.4e1 - t662 * t790 * t912 + (-t924 - t927) * (-t613 - (-t876 / 0.2e1 - t873 / 0.2e1 - t813 / 0.4e1) * t689) - (t673 * t906 + t921) * t597) * t662) * t869;
t663 = t662 * t904;
t771 = t813 + 0.2e1 * t873 + 0.2e1 * t876;
t798 = -t862 / 0.2e1;
t588 = (-t663 * t910 - (-t597 + t662 * t615 * t921 + (t613 * t662 + t771 * t798) * t673) * t615) * t689;
t850 = t725 * t729;
t603 = t615 * t850;
t651 = t723 * mrSges(3,2) - Ifges(3,6);
t654 = t723 * mrSges(3,1) - Ifges(3,5);
t618 = -(-t651 * t735 - t654 * t729) * t718 - t717 * (t651 * t729 - t654 * t735);
t889 = mrSges(3,2) * t729;
t778 = mrSges(3,1) * t735 - t889;
t683 = t729 * mrSges(3,1);
t686 = mrSges(3,2) * t735;
t838 = t686 + t683;
t621 = -(-t778 - t881) * t718 - (mrSges(2,2) + t838) * t717 + t695;
t639 = t742 * qJ(2,2) + t818;
t714 = t735 ^ 2;
t853 = t725 * t714;
t868 = t615 * t717;
t783 = t853 * t868;
t786 = t663 * t667 * t855;
t808 = t662 * t845;
t823 = pkin(2) * t889;
t816 = (t823 - t725) * t868;
t826 = mrSges(3,1) * t613;
t886 = Ifges(3,4) * t729;
t832 = 0.4e1 * t886;
t847 = t729 * t735;
t883 = t714 * Ifges(3,4);
t915 = 0.2e1 * t883 - Ifges(3,4);
t841 = -(-((t710 + t832) * t735 - 0.2e1 * t823 + t780) * t702 - (t714 * t834 + t791 * t735 + (t881 - t889) * t906) * t718 - t853 + (t683 + mrSges(2,2)) * t835 - t723 ^ 2 * m(3) - m(2) * qJ(2,2) ^ 2 + t776) * t588 + t621 * t582 - t618 * t786 - 0.2e1 * (Ifges(3,4) * t847 + t702 * t853 - t726 * qJ(2,2) + (-(-pkin(2) * t683 + t725 * t847 + t893) * t718 + pkin(1) * t686) * t717) * t588 - 0.2e1 * ((((-qJ(2,2) / 0.2e1 + t755) * mrSges(3,2) + t743) * t808 - t826) * t717 - t603) * t735 * t808 - 0.2e1 * t615 * (Ifges(3,4) * t808 - t639 * t597) - 0.4e1 * (((t693 + t850) * t735 + t729 * t820 + t915) * t702 - t883) * t615 * t808 + ((-t783 + (((-qJ(2,2) / 0.4e1 + t754) * mrSges(3,1) + t744) * t808 - ((t694 + t886) * t901 + t821) * t615) * t735 - t816 / 0.2e1 - t729 * (-t651 * t808 - 0.2e1 * t826) / 0.4e1) * t795 - (-t615 * t836 - t654 * t808) * t729 * t805) * t662;
t724 = pkin(5) + qJ(2,1);
t698 = -pkin(6) - t724;
t690 = 0.1e1 / t698;
t646 = sin(t659);
t668 = sin(t708);
t731 = sin(qJ(3,1));
t625 = t668 * t906 + pkin(3) * t646 + (sin(t705) + t731) * pkin(2);
t674 = cos(t708);
t664 = 0.1e1 / t674;
t800 = t625 * t664 / 0.2e1;
t732 = sin(qJ(1,1));
t738 = cos(qJ(1,1));
t628 = t732 * t656 + t738 * t698;
t631 = t656 * t738 - t732 * t698;
t721 = legFrame(1,3);
t678 = sin(t721);
t681 = cos(t721);
t637 = t678 * t738 + t681 * t732;
t894 = pkin(3) * t674;
t610 = t628 * t681 + t631 * t678 + t637 * t894;
t872 = t610 * t740;
t634 = -t678 * t732 + t681 * t738;
t607 = -t628 * t678 + t631 * t681 + t634 * t894;
t875 = t607 * t741;
t598 = (t739 * t800 + t872 + t875) * t690;
t858 = t664 * t739;
t616 = (t634 * t741 + t637 * t740 + t668 * t858) * t690;
t611 = t616 * pkin(1);
t758 = 2 * qJ(3,1);
t812 = t625 * t858;
t859 = t664 * t690;
t867 = t616 * t690;
t903 = 0.1e1 / t674 ^ 2;
t911 = t646 * t698;
t922 = t928 / 0.2e1 + t925 / 0.2e1;
t583 = (t656 + t894) * t903 * t910 * t859 + (t903 * t897 * t911 + (-(-t759 * cos((3 * t708)) + (-0.4e1 * t698 ^ 2 + t803) * t674 + t804 + (cos((t752 + t758)) + cos((pkin(7) + t758))) * t843 + (-cos((-pkin(7) + qJ(3,1))) - cos((t752 + qJ(3,1)))) * t762) * t616 / 0.4e1 - t664 * t790 * t911 + (-t925 - t928) * (-t611 - (-t875 / 0.2e1 - t872 / 0.2e1 - t812 / 0.4e1) * t690) - (t674 * t906 + t922) * t598) * t664) * t867;
t665 = t664 * t903;
t770 = t812 + 0.2e1 * t872 + 0.2e1 * t875;
t797 = -t859 / 0.2e1;
t589 = (-t665 * t910 - (-t598 + t664 * t616 * t922 + (t611 * t664 + t770 * t797) * t674) * t616) * t690;
t849 = t725 * t731;
t604 = t616 * t849;
t652 = t724 * mrSges(3,2) - Ifges(3,6);
t655 = t724 * mrSges(3,1) - Ifges(3,5);
t619 = -(-t652 * t737 - t655 * t731) * t718 - t717 * (t652 * t731 - t655 * t737);
t888 = mrSges(3,2) * t731;
t777 = mrSges(3,1) * t737 - t888;
t684 = t731 * mrSges(3,1);
t687 = mrSges(3,2) * t737;
t837 = t687 + t684;
t622 = -(-t777 - t881) * t718 - (mrSges(2,2) + t837) * t717 + t695;
t640 = t742 * qJ(2,1) + t818;
t715 = t737 ^ 2;
t852 = t725 * t715;
t866 = t616 * t717;
t782 = t852 * t866;
t785 = t665 * t668 * t855;
t806 = t664 * t845;
t822 = pkin(2) * t888;
t815 = (t822 - t725) * t866;
t825 = mrSges(3,1) * t611;
t885 = Ifges(3,4) * t731;
t831 = 0.4e1 * t885;
t846 = t731 * t737;
t882 = t715 * Ifges(3,4);
t916 = 0.2e1 * t882 - Ifges(3,4);
t840 = -(-((t710 + t831) * t737 - 0.2e1 * t822 + t780) * t702 - (t715 * t834 + t791 * t737 + (t881 - t888) * t906) * t718 - t852 + (t684 + mrSges(2,2)) * t835 - t724 ^ 2 * m(3) - m(2) * qJ(2,1) ^ 2 + t776) * t589 + t622 * t583 - t619 * t785 - 0.2e1 * (Ifges(3,4) * t846 + t702 * t852 - t726 * qJ(2,1) + (-(-pkin(2) * t684 + t725 * t846 + t893) * t718 + pkin(1) * t687) * t717) * t589 - 0.2e1 * ((((-qJ(2,1) / 0.2e1 + t755) * mrSges(3,2) + t743) * t806 - t825) * t717 - t604) * t737 * t806 - 0.2e1 * t616 * (Ifges(3,4) * t806 - t640 * t598) - 0.4e1 * (((t693 + t849) * t737 + t731 * t820 + t916) * t702 - t882) * t616 * t806 + ((-t782 + (((-qJ(2,1) / 0.4e1 + t754) * mrSges(3,1) + t744) * t806 - ((t694 + t885) * t901 + t821) * t616) * t737 - t815 / 0.2e1 - t731 * (-t652 * t806 - 0.2e1 * t825) / 0.4e1) * t795 - (-t616 * t836 - t655 * t806) * t731 * t805) * t664;
t902 = -0.4e1 * t702;
t892 = mrSges(3,1) * t717;
t891 = mrSges(3,2) * t717;
t593 = t770 * t690 - t611;
t880 = t593 * t731;
t594 = t772 * t688 - t612;
t879 = t594 * t727;
t595 = t771 * t689 - t613;
t878 = t595 * t729;
t819 = -t899 / 0.4e1;
t796 = 0.2e1 * t845;
t578 = -t742 * t581 - t620 * t587;
t599 = -t638 * t614 + (-t779 * t717 - t839 * t718) * t660 * t796;
t794 = t599 * t614 + t578;
t579 = -t742 * t582 - t621 * t588;
t600 = -t639 * t615 + (-t778 * t717 - t838 * t718) * t662 * t796;
t793 = t600 * t615 + t579;
t580 = -t742 * t583 - t622 * t589;
t601 = -t640 * t616 + (-t777 * t717 - t837 * t718) * t664 * t796;
t792 = t601 * t616 + t580;
t745 = -Ifges(3,4) / 0.2e1;
t692 = -t898 / 0.4e1;
t691 = Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1;
t1 = [-(t792 * t607 + t840 * t634) * t690 - (t793 * t606 + t841 * t633) * t689 - (t794 * t605 + t842 * t632) * t688; -(t792 * t610 + t840 * t637) * t690 - (t793 * t609 + t841 * t636) * t689 - (t794 * t608 + t842 * t635) * t688; t623 * t578 * t799 + t624 * t579 * t798 + t625 * t580 * t797 - t599 * t802 * t871 - t600 * t801 * t869 - t601 * t800 * t867 - t842 * t666 * t865 - t841 * t667 * t862 - t840 * t668 * t859 + ((Ifges(3,3) * t787 - t617 * t587 - t614 * ((0.2e1 * t784 + (mrSges(3,2) * t594 - (t833 + t899) * t870) * t733 + t817 + mrSges(3,1) * t879) * t718 + (t594 * t892 - t602) * t733 - t879 * t891 - ((t884 + (t691 * t727 + t692) * t733 + t727 * t819 + t745) * t902 + t914) * t614)) * t660 + (Ifges(3,3) * t786 - t618 * t588 - t615 * ((0.2e1 * t783 + (mrSges(3,2) * t595 - (t832 + t899) * t868) * t735 + t816 + mrSges(3,1) * t878) * t718 + (t595 * t892 - t603) * t735 - t878 * t891 - ((t883 + (t691 * t729 + t692) * t735 + t729 * t819 + t745) * t902 + t915) * t615)) * t662 + (Ifges(3,3) * t785 - t619 * t589 - t616 * ((0.2e1 * t782 + (mrSges(3,2) * t593 - (t831 + t899) * t866) * t737 + t815 + mrSges(3,1) * t880) * t718 + (t593 * t892 - t604) * t737 - t880 * t891 - ((t882 + (t691 * t731 + t692) * t737 + t731 * t819 + t745) * t902 + t916) * t616)) * t664) * t760;];
taucX  = t1;
