% Calculate vector of centrifugal and Coriolis load for parallel robot
% P3RRPRR8V2G1A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:11:15
% EndTime: 2022-11-07 13:11:21
% DurationCPUTime: 5.82s
% Computational Cost: add. (26127->431), mult. (28620->667), div. (4992->9), fcn. (20493->62), ass. (0->304)
t781 = qJ(3,1) + pkin(5);
t697 = mrSges(3,2) * t781 - Ifges(3,6);
t774 = sin(pkin(7));
t775 = cos(pkin(7));
t957 = mrSges(3,1) * t781 - Ifges(3,5);
t967 = -t697 * t774 + t957 * t775;
t780 = qJ(3,2) + pkin(5);
t696 = mrSges(3,2) * t780 - Ifges(3,6);
t958 = mrSges(3,1) * t780 - Ifges(3,5);
t966 = -t696 * t774 + t958 * t775;
t779 = qJ(3,3) + pkin(5);
t695 = mrSges(3,2) * t779 - Ifges(3,6);
t959 = mrSges(3,1) * t779 - Ifges(3,5);
t965 = -t695 * t774 + t959 * t775;
t782 = Ifges(3,2) - Ifges(3,1);
t686 = -pkin(2) * mrSges(3,2) - t774 * t782;
t927 = mrSges(3,1) * t774;
t688 = -pkin(2) * t927 + Ifges(2,4) - Ifges(3,4);
t756 = t775 ^ 2;
t925 = Ifges(3,4) * t756;
t656 = t686 * t775 + t688 + 0.2e1 * t925;
t798 = xDP(3);
t829 = t656 * t798;
t964 = 2 * pkin(1);
t963 = -0.2e1 * pkin(2);
t752 = -pkin(6) - t781;
t741 = 0.1e1 / t752;
t762 = qJ(2,1) + pkin(7);
t789 = sin(qJ(2,1));
t684 = pkin(2) * t789 + pkin(3) * sin(t762);
t814 = 0.2e1 * qJ(2,1);
t761 = t814 + pkin(7);
t718 = sin(t761);
t944 = pkin(2) * pkin(3);
t892 = 0.2e1 * t944;
t819 = pkin(2) ^ 2;
t895 = t819 * sin(t814);
t707 = 0.2e1 * t762;
t817 = pkin(3) ^ 2;
t898 = t817 * sin(t707);
t835 = t718 * t892 + t895 + t898;
t651 = t684 * t964 + t835;
t795 = cos(qJ(2,1));
t748 = t795 * pkin(2);
t730 = cos(t762);
t930 = pkin(3) * t730;
t681 = t748 + t930;
t676 = 0.1e1 / t681;
t904 = t676 * t798;
t870 = t651 * t904;
t710 = t748 + pkin(1);
t790 = sin(qJ(1,1));
t796 = cos(qJ(1,1));
t659 = t710 * t790 + t752 * t796;
t662 = t710 * t796 - t752 * t790;
t778 = legFrame(1,3);
t735 = sin(t778);
t738 = cos(t778);
t669 = -t735 * t790 + t738 * t796;
t645 = -t659 * t735 + t662 * t738 + t669 * t930;
t800 = xDP(1);
t913 = t645 * t800;
t668 = t735 * t796 + t738 * t790;
t644 = t659 * t738 + t662 * t735 + t668 * t930;
t799 = xDP(2);
t914 = t644 * t799;
t621 = (t913 + t914 + t870 / 0.2e1) * t741;
t751 = -pkin(6) - t780;
t740 = 0.1e1 / t751;
t760 = qJ(2,2) + pkin(7);
t787 = sin(qJ(2,2));
t683 = pkin(2) * t787 + pkin(3) * sin(t760);
t812 = 0.2e1 * qJ(2,2);
t759 = t812 + pkin(7);
t716 = sin(t759);
t896 = t819 * sin(t812);
t706 = 0.2e1 * t760;
t899 = t817 * sin(t706);
t836 = t716 * t892 + t896 + t899;
t650 = t683 * t964 + t836;
t793 = cos(qJ(2,2));
t747 = t793 * pkin(2);
t726 = cos(t760);
t931 = pkin(3) * t726;
t680 = t747 + t931;
t674 = 0.1e1 / t680;
t906 = t674 * t798;
t871 = t650 * t906;
t709 = t747 + pkin(1);
t788 = sin(qJ(1,2));
t794 = cos(qJ(1,2));
t658 = t709 * t788 + t751 * t794;
t661 = t709 * t794 - t751 * t788;
t777 = legFrame(2,3);
t734 = sin(t777);
t737 = cos(t777);
t667 = -t734 * t788 + t737 * t794;
t643 = -t658 * t734 + t661 * t737 + t667 * t931;
t915 = t643 * t800;
t666 = t734 * t794 + t737 * t788;
t642 = t658 * t737 + t661 * t734 + t666 * t931;
t916 = t642 * t799;
t620 = (t915 + t916 + t871 / 0.2e1) * t740;
t750 = -pkin(6) - t779;
t739 = 0.1e1 / t750;
t758 = qJ(2,3) + pkin(7);
t785 = sin(qJ(2,3));
t682 = pkin(2) * t785 + pkin(3) * sin(t758);
t810 = 0.2e1 * qJ(2,3);
t757 = t810 + pkin(7);
t714 = sin(t757);
t897 = t819 * sin(t810);
t705 = 0.2e1 * t758;
t900 = t817 * sin(t705);
t837 = t714 * t892 + t897 + t900;
t649 = t682 * t964 + t837;
t791 = cos(qJ(2,3));
t746 = t791 * pkin(2);
t722 = cos(t758);
t932 = pkin(3) * t722;
t679 = t746 + t932;
t672 = 0.1e1 / t679;
t908 = t672 * t798;
t872 = t649 * t908;
t708 = t746 + pkin(1);
t786 = sin(qJ(1,3));
t792 = cos(qJ(1,3));
t657 = t708 * t786 + t750 * t792;
t660 = t708 * t792 - t750 * t786;
t776 = legFrame(3,3);
t733 = sin(t776);
t736 = cos(t776);
t665 = -t733 * t786 + t736 * t792;
t641 = -t657 * t733 + t660 * t736 + t665 * t932;
t917 = t641 * t800;
t664 = t733 * t792 + t736 * t786;
t640 = t657 * t736 + t660 * t733 + t664 * t932;
t918 = t640 * t799;
t619 = (t917 + t918 + t872 / 0.2e1) * t739;
t634 = (t664 * t799 + t665 * t800 + t682 * t908) * t739;
t632 = pkin(1) * t634;
t611 = -t632 - (-t917 / 0.2e1 - t918 / 0.2e1 - t872 / 0.4e1) * t739;
t820 = pkin(1) ^ 2;
t625 = (t750 ^ 2 + t820) * t634;
t646 = t750 * t908;
t807 = 0.2e1 * pkin(7);
t720 = cos(t807 + qJ(2,3));
t723 = cos(qJ(2,3) - pkin(7));
t773 = t798 ^ 2;
t809 = 0.3e1 * qJ(2,3);
t816 = pkin(3) * t817;
t818 = pkin(2) * t819;
t749 = t817 + t819;
t844 = t817 * cos(t705) + t819 * cos(t810) + t749;
t850 = 0.3e1 / 0.4e1 * t819;
t851 = 0.3e1 / 0.4e1 * t817;
t929 = pkin(3) * t819;
t855 = -0.2e1 * t816 - 0.4e1 * t929;
t909 = t672 * t739;
t861 = -t909 / 0.2e1;
t880 = -0.2e1 * t929;
t933 = pkin(2) * t817;
t881 = -0.2e1 * t933;
t885 = -0.6e1 * t819 - 0.3e1 * t817;
t721 = cos(t757);
t888 = t721 + t775;
t936 = -0.3e1 / 0.4e1 * t819;
t937 = -t798 / 0.2e1;
t940 = -t634 / 0.4e1;
t945 = -0.2e1 * t818 - 0.4e1 * t933;
t879 = t775 * t944;
t946 = -0.4e1 * pkin(1) * (t879 + t819 / 0.2e1 + t817 / 0.2e1);
t953 = 0.1e1 / t679 ^ 2;
t601 = (t720 * t881 + t855 * t722 + t723 * t880 + t791 * t945 + t946) * t773 * t953 * t861 + (t837 * t953 * t750 * t937 + ((-t816 * cos(0.3e1 * t758) - t818 * cos(t809)) * t940 - (t900 / 0.2e1 + t897 / 0.2e1) * t646 + (-(-cos(pkin(7) + t809) - t723) * t634 * t850 + (-t619 * t964 + t885 * t940 + t625) * t722) * pkin(3) + ((-t634 * t936 + t625) * t791 - (-cos(t807 + t809) - t720 - 0.2e1 * t791) * t634 * t851 + (-0.2e1 * t611 * t888 - t646 * t714) * pkin(3) - (t888 * pkin(3) + t791 * t964) * t619) * pkin(2) + (-t619 / 0.2e1 - t611) * t844) * t672) * t739 * t634;
t673 = t672 * t953;
t863 = -0.2e1 * t879;
t875 = t634 * t909;
t893 = -0.2e1 * t944;
t903 = (0.2e1 * t879 + t749) * t773;
t604 = t673 * t739 * t903 + (-(t721 * t893 - t844 + t863) * t634 / 0.2e1 + (-0.2e1 * t619 + t632) * t679) * t875;
t805 = m(3) * pkin(2);
t742 = mrSges(2,1) + t805;
t934 = pkin(1) * t742;
t628 = t634 * t934;
t856 = -m(3) * qJ(3,3) - mrSges(3,3);
t941 = pkin(5) * mrSges(2,2);
t862 = Ifges(2,6) - t941;
t922 = -t742 * pkin(5) + Ifges(2,5);
t637 = -(t856 * pkin(2) + t922 - t965) * t785 - t791 * (-t695 * t775 - t774 * t959 + t862);
t926 = mrSges(3,2) * t774;
t846 = mrSges(3,1) * t775 - t926;
t670 = -t805 - t846;
t806 = m(3) * pkin(1);
t845 = mrSges(3,2) * t775 + t927;
t652 = -t670 * t791 - t785 * t845 + t806;
t901 = t782 * t756;
t924 = Ifges(3,4) * t774;
t935 = m(3) * t819;
t955 = 0.2e1 * mrSges(3,1);
t655 = 0.2e1 * t901 + (pkin(2) * t955 + 0.4e1 * t924) * t775 + t935 + t926 * t963 + Ifges(2,2) - Ifges(2,1) - t782;
t678 = 0.2e1 * t686;
t685 = (t742 - t926) * t964;
t804 = m(3) * pkin(5);
t701 = t804 - t856;
t808 = -pkin(5) / 0.2e1;
t743 = -qJ(3,3) / 0.2e1 + t808;
t769 = pkin(1) * t955;
t770 = t791 ^ 2;
t802 = Ifges(3,6) / 0.2e1;
t803 = Ifges(3,5) / 0.2e1;
t830 = -m(2) * pkin(5) ^ 2 - (m(2) + m(3)) * t820 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t901;
t834 = (mrSges(2,2) + t845) * t964;
t868 = t673 * t682 * t773;
t876 = 0.4e1 * t829;
t878 = 0.4e1 * t925;
t882 = (mrSges(2,2) + t927) * t964;
t902 = (-t941 / 0.2e1 + Ifges(2,6) / 0.2e1) * t798;
t912 = t655 * t785;
t921 = t634 * t672;
t923 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t928 = (mrSges(2,3) + mrSges(3,3)) * pkin(5);
t942 = mrSges(3,2) * pkin(1);
t950 = -0.2e1 * t634;
t954 = -0.2e1 * mrSges(3,3);
t891 = (-t655 * t770 - (t785 * t878 + (t678 * t785 + t769) * t775 + t685) * t791 + t785 * t882 - t779 ^ 2 * m(3) + qJ(3,3) * t954 + t830) * t604 - t637 * t868 + t652 * t601 + 0.2e1 * (-t688 * t785 * t791 - (-t785 * t942 - t924) * t775 - t928) * t604 - t770 * t876 * t921 - ((t701 * pkin(2) + t923 + t965) * t908 - (t834 + 0.2e1 * t912) * t634) * t791 * t908 - 0.2e1 * (((mrSges(3,2) * t743 + t802) * t908 - mrSges(3,1) * t632) * t775 + ((mrSges(3,1) * t743 + t803) * t908 + mrSges(3,2) * t632) * t774 + t672 * t902 - t628) * t785 * t908 + (-t619 * t701 - t672 * t829) * t950;
t635 = (t666 * t799 + t667 * t800 + t683 * t906) * t740;
t633 = pkin(1) * t635;
t612 = -t633 - (-t915 / 0.2e1 - t916 / 0.2e1 - t871 / 0.4e1) * t740;
t626 = (t751 ^ 2 + t820) * t635;
t647 = t751 * t906;
t725 = cos(qJ(2,2) + t807);
t727 = cos(qJ(2,2) - pkin(7));
t811 = 0.3e1 * qJ(2,2);
t843 = t817 * cos(t706) + t819 * cos(t812) + t749;
t907 = t674 * t740;
t860 = -t907 / 0.2e1;
t724 = cos(t759);
t887 = t724 + t775;
t939 = -t635 / 0.4e1;
t952 = 0.1e1 / t680 ^ 2;
t602 = (t725 * t881 + t855 * t726 + t727 * t880 + t793 * t945 + t946) * t773 * t952 * t860 + (t836 * t952 * t751 * t937 + ((-t816 * cos(0.3e1 * t760) - t818 * cos(t811)) * t939 - (t899 / 0.2e1 + t896 / 0.2e1) * t647 + (-(-cos(t811 + pkin(7)) - t727) * t635 * t850 + (-t620 * t964 + t885 * t939 + t626) * t726) * pkin(3) + ((-t635 * t936 + t626) * t793 - (-cos(t811 + t807) - t725 - 0.2e1 * t793) * t635 * t851 + (-0.2e1 * t612 * t887 - t647 * t716) * pkin(3) - (t887 * pkin(3) + t793 * t964) * t620) * pkin(2) + (-t620 / 0.2e1 - t612) * t843) * t674) * t740 * t635;
t675 = t674 * t952;
t874 = t635 * t907;
t605 = t675 * t740 * t903 + (-(t724 * t893 - t843 + t863) * t635 / 0.2e1 + (-0.2e1 * t620 + t633) * t680) * t874;
t629 = t635 * t934;
t857 = -m(3) * qJ(3,2) - mrSges(3,3);
t638 = -(t857 * pkin(2) + t922 - t966) * t787 - t793 * (-t696 * t775 - t774 * t958 + t862);
t653 = -t670 * t793 - t787 * t845 + t806;
t702 = t804 - t857;
t744 = -qJ(3,2) / 0.2e1 + t808;
t771 = t793 ^ 2;
t866 = t675 * t683 * t773;
t911 = t655 * t787;
t920 = t635 * t674;
t949 = -0.2e1 * t635;
t890 = (-t655 * t771 - (t787 * t878 + (t678 * t787 + t769) * t775 + t685) * t793 + t787 * t882 - t780 ^ 2 * m(3) + qJ(3,2) * t954 + t830) * t605 - t638 * t866 + t653 * t602 + 0.2e1 * (-t688 * t787 * t793 - (-t787 * t942 - t924) * t775 - t928) * t605 - t771 * t876 * t920 - ((t702 * pkin(2) + t923 + t966) * t906 - (t834 + 0.2e1 * t911) * t635) * t793 * t906 - 0.2e1 * (((mrSges(3,2) * t744 + t802) * t906 - mrSges(3,1) * t633) * t775 + ((mrSges(3,1) * t744 + t803) * t906 + mrSges(3,2) * t633) * t774 + t674 * t902 - t629) * t787 * t906 + (-t620 * t702 - t674 * t829) * t949;
t636 = (t668 * t799 + t669 * t800 + t684 * t904) * t741;
t631 = t636 * pkin(1);
t610 = -t631 - (-t913 / 0.2e1 - t914 / 0.2e1 - t870 / 0.4e1) * t741;
t627 = (t752 ^ 2 + t820) * t636;
t648 = t752 * t904;
t729 = cos(qJ(2,1) + t807);
t731 = cos(qJ(2,1) - pkin(7));
t813 = 0.3e1 * qJ(2,1);
t842 = cos(t707) * t817 + cos(t814) * t819 + t749;
t905 = t676 * t741;
t859 = -t905 / 0.2e1;
t728 = cos(t761);
t886 = t728 + t775;
t938 = -t636 / 0.4e1;
t951 = 0.1e1 / t681 ^ 2;
t603 = (t729 * t881 + t855 * t730 + t731 * t880 + t795 * t945 + t946) * t773 * t951 * t859 + (t835 * t951 * t752 * t937 + ((-t816 * cos(0.3e1 * t762) - t818 * cos(t813)) * t938 - (t898 / 0.2e1 + t895 / 0.2e1) * t648 + (-(-cos(pkin(7) + t813) - t731) * t636 * t850 + (-t621 * t964 + t885 * t938 + t627) * t730) * pkin(3) + ((-t636 * t936 + t627) * t795 - (-cos(t813 + t807) - t729 - 0.2e1 * t795) * t636 * t851 + (-0.2e1 * t610 * t886 - t648 * t718) * pkin(3) - (t886 * pkin(3) + t795 * t964) * t621) * pkin(2) + (-t621 / 0.2e1 - t610) * t842) * t676) * t741 * t636;
t677 = t676 * t951;
t873 = t636 * t905;
t606 = t677 * t741 * t903 + (-(t728 * t893 - t842 + t863) * t636 / 0.2e1 + (-0.2e1 * t621 + t631) * t681) * t873;
t630 = t742 * t631;
t858 = -m(3) * qJ(3,1) - mrSges(3,3);
t639 = -(t858 * pkin(2) + t922 - t967) * t789 - t795 * (-t697 * t775 - t774 * t957 + t862);
t654 = -t670 * t795 - t789 * t845 + t806;
t703 = t804 - t858;
t745 = -qJ(3,1) / 0.2e1 + t808;
t772 = t795 ^ 2;
t864 = t677 * t684 * t773;
t910 = t655 * t789;
t919 = t636 * t676;
t948 = -0.2e1 * t636;
t889 = (-t655 * t772 - (t789 * t878 + (t678 * t789 + t769) * t775 + t685) * t795 + t789 * t882 - t781 ^ 2 * m(3) + qJ(3,1) * t954 + t830) * t606 - t639 * t864 + t654 * t603 + 0.2e1 * (-t688 * t789 * t795 - (-t789 * t942 - t924) * t775 - t928) * t606 - t772 * t876 * t919 - ((t703 * pkin(2) + t923 + t967) * t904 - (t834 + 0.2e1 * t910) * t636) * t795 * t904 - 0.2e1 * (((mrSges(3,2) * t745 + t802) * t904 - mrSges(3,1) * t631) * t775 + ((mrSges(3,1) * t745 + t803) * t904 + mrSges(3,2) * t631) * t774 + t676 * t902 - t630) * t789 * t904 + (-t621 * t703 - t676 * t829) * t948;
t947 = -0.2e1 * t656;
t943 = mrSges(2,2) * pkin(1);
t894 = -0.2e1 * t805;
t598 = -m(3) * t601 + t604 * t652;
t841 = t845 * t791;
t622 = (m(3) * t779 + mrSges(3,3)) * t634 / 0.2e1 + (-t670 * t785 + t841) * t908;
t854 = t622 * t950 + t598;
t599 = -m(3) * t602 + t605 * t653;
t840 = t845 * t793;
t623 = (m(3) * t780 + mrSges(3,3)) * t635 / 0.2e1 + (-t670 * t787 + t840) * t906;
t853 = t623 * t949 + t599;
t600 = -m(3) * t603 + t606 * t654;
t839 = t845 * t795;
t624 = (m(3) * t781 + mrSges(3,3)) * t636 / 0.2e1 + (-t670 * t789 + t839) * t904;
t852 = t624 * t948 + t600;
t663 = t846 * t963 - Ifges(2,3) - Ifges(3,3) - t935;
t1 = [-(t852 * t645 + t889 * t669) * t741 - (t853 * t643 + t890 * t667) * t740 - (t854 * t641 + t891 * t665) * t739; -(t852 * t644 + t889 * t668) * t741 - (t853 * t642 + t890 * t666) * t740 - (t854 * t640 + t891 * t664) * t739; t676 * (t606 * t639 - t663 * t864) - ((-t621 * t894 - t630) * t789 + (t789 * t846 + t839) * (-t631 - (-t870 - 0.2e1 * t913 - 0.2e1 * t914) * t741) - (t772 * t947 + (t910 + t943) * t795 + t656) * t636) * t919 + t674 * (t605 * t638 - t663 * t866) - ((-t620 * t894 - t629) * t787 + (t787 * t846 + t840) * (-t633 - (-t871 - 0.2e1 * t915 - 0.2e1 * t916) * t740) - (t771 * t947 + (t911 + t943) * t793 + t656) * t635) * t920 + t672 * (t604 * t637 - t663 * t868) - ((-t619 * t894 - t628) * t785 + (t785 * t846 + t841) * (-t632 - (-t872 - 0.2e1 * t917 - 0.2e1 * t918) * t739) - (t770 * t947 + (t912 + t943) * t791 + t656) * t634) * t921 - t891 * t682 * t909 - t890 * t683 * t907 - t889 * t684 * t905 + (t600 * t859 + t624 * t873) * t651 + (t599 * t860 + t623 * t874) * t650 + (t598 * t861 + t622 * t875) * t649;];
taucX  = t1;
