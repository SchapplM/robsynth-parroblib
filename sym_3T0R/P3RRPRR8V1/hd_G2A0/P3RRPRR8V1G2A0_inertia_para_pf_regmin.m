% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x15]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:48
% EndTime: 2022-11-04 17:04:54
% DurationCPUTime: 5.76s
% Computational Cost: add. (6763->439), mult. (12885->995), div. (1527->18), fcn. (11964->23), ass. (0->423)
t1058 = 2 * qJ(3,1);
t1057 = 2 * qJ(3,2);
t1056 = 2 * qJ(3,3);
t790 = cos(pkin(5));
t764 = pkin(2) * t790 + pkin(1);
t798 = sin(qJ(1,3));
t1013 = t764 * t798;
t789 = sin(pkin(5));
t1047 = pkin(2) * t789;
t794 = legFrame(3,2);
t768 = sin(t794);
t771 = cos(t794);
t797 = sin(qJ(2,3));
t803 = cos(qJ(2,3));
t973 = t798 * t1047;
t713 = (-t768 * t1013 + t771 * t1047) * t803 + (t764 * t771 + t768 * t973) * t797;
t716 = (t771 * t1013 + t768 * t1047) * t803 + (t764 * t768 - t771 * t973) * t797;
t1055 = t713 * t716;
t800 = sin(qJ(1,2));
t1012 = t764 * t800;
t795 = legFrame(2,2);
t769 = sin(t795);
t772 = cos(t795);
t799 = sin(qJ(2,2));
t805 = cos(qJ(2,2));
t972 = t800 * t1047;
t714 = (-t769 * t1012 + t772 * t1047) * t805 + (t764 * t772 + t769 * t972) * t799;
t717 = (t772 * t1012 + t769 * t1047) * t805 + (t764 * t769 - t772 * t972) * t799;
t1054 = t714 * t717;
t802 = sin(qJ(1,1));
t1011 = t764 * t802;
t796 = legFrame(1,2);
t770 = sin(t796);
t773 = cos(t796);
t801 = sin(qJ(2,1));
t807 = cos(qJ(2,1));
t971 = t802 * t1047;
t715 = (-t770 * t1011 + t773 * t1047) * t807 + (t764 * t773 + t770 * t971) * t801;
t718 = (t773 * t1011 + t770 * t1047) * t807 + (t764 * t770 - t773 * t971) * t801;
t1053 = t715 * t718;
t1052 = -0.2e1 * pkin(1);
t1051 = pkin(1) * t790;
t1050 = pkin(1) * t803;
t1049 = pkin(1) * t805;
t1048 = pkin(1) * t807;
t995 = t789 * t797;
t746 = -pkin(2) * t995 + t764 * t803;
t737 = 0.1e1 / t746;
t1046 = t713 * t737;
t791 = pkin(4) + qJ(3,3);
t774 = 1 / t791;
t1045 = t713 * t774;
t994 = t789 * t799;
t747 = -pkin(2) * t994 + t764 * t805;
t739 = 0.1e1 / t747;
t1044 = t714 * t739;
t792 = pkin(4) + qJ(3,2);
t776 = 1 / t792;
t1043 = t714 * t776;
t993 = t789 * t801;
t748 = -pkin(2) * t993 + t764 * t807;
t741 = 0.1e1 / t748;
t1042 = t715 * t741;
t793 = pkin(4) + qJ(3,1);
t778 = 1 / t793;
t1041 = t715 * t778;
t1040 = t716 * t737;
t1039 = t716 * t774;
t1038 = t717 * t739;
t1037 = t717 * t776;
t1036 = t718 * t741;
t1035 = t718 * t778;
t761 = pkin(2) * cos(qJ(2,3) + pkin(5)) + t1050;
t804 = cos(qJ(1,3));
t734 = t761 * t804 + t791 * t798;
t783 = t803 ^ 2;
t812 = pkin(1) ^ 2;
t907 = (qJ(3,3) ^ 2) + t783 * t812;
t719 = (-t734 * t1050 + t907 * t804) * t774;
t1034 = t719 * t737;
t762 = pkin(2) * cos(qJ(2,2) + pkin(5)) + t1049;
t806 = cos(qJ(1,2));
t735 = t762 * t806 + t792 * t800;
t785 = t805 ^ 2;
t906 = (qJ(3,2) ^ 2) + t785 * t812;
t720 = (-t735 * t1049 + t906 * t806) * t776;
t1033 = t720 * t739;
t763 = pkin(2) * cos(qJ(2,1) + pkin(5)) + t1048;
t808 = cos(qJ(1,1));
t736 = t763 * t808 + t793 * t802;
t787 = t807 ^ 2;
t905 = (qJ(3,1) ^ 2) + t787 * t812;
t721 = (-t736 * t1048 + t905 * t808) * t778;
t1032 = t721 * t741;
t1031 = t737 * t774;
t738 = 0.1e1 / t746 ^ 2;
t775 = 1 / t791 ^ 2;
t1030 = t738 * t775;
t1029 = t739 * t776;
t740 = 0.1e1 / t747 ^ 2;
t777 = 1 / t792 ^ 2;
t1028 = t740 * t777;
t1027 = t741 * t778;
t742 = 0.1e1 / t748 ^ 2;
t779 = 1 / t793 ^ 2;
t1026 = t742 * t779;
t749 = -t790 * t803 + t995;
t1025 = t749 * t775;
t750 = -t790 * t805 + t994;
t1024 = t750 * t777;
t751 = -t790 * t807 + t993;
t1023 = t751 * t779;
t989 = t790 * t797;
t992 = t789 * t803;
t752 = t989 + t992;
t1022 = t752 * t775;
t988 = t790 * t799;
t991 = t789 * t805;
t753 = t988 + t991;
t1021 = t753 * t777;
t987 = t790 * t801;
t990 = t789 * t807;
t754 = t987 + t990;
t1020 = t754 * t779;
t755 = 0.1e1 / t761;
t1019 = t755 * t768;
t1018 = t755 * t771;
t757 = 0.1e1 / t762;
t1017 = t757 * t769;
t1016 = t757 * t772;
t759 = 0.1e1 / t763;
t1015 = t759 * t770;
t1014 = t759 * t773;
t1010 = t774 * t789;
t1009 = t774 * t790;
t1008 = t774 * t804;
t784 = t804 ^ 2;
t1007 = t775 * t784;
t1006 = t775 * t804;
t1005 = t776 * t789;
t1004 = t776 * t790;
t1003 = t776 * t806;
t786 = t806 ^ 2;
t1002 = t777 * t786;
t1001 = t777 * t806;
t1000 = t778 * t789;
t999 = t778 * t790;
t998 = t778 * t808;
t788 = t808 ^ 2;
t997 = t779 * t788;
t996 = t779 * t808;
t986 = t797 * t803;
t985 = t799 * t805;
t984 = t801 * t807;
t983 = 0.2e1 * pkin(1);
t982 = t789 * t1052;
t981 = -0.2e1 * t1051;
t980 = 0.2e1 * t1051;
t979 = t737 * t1050;
t978 = t739 * t1049;
t977 = t741 * t1048;
t976 = pkin(1) * t755 * t797;
t975 = pkin(1) * t757 * t799;
t974 = pkin(1) * t759 * t801;
t970 = qJ(3,1) * t759 * t789;
t969 = qJ(3,2) * t757 * t789;
t968 = qJ(3,3) * t755 * t789;
t725 = t746 * t798 - t791 * t804;
t743 = pkin(2) * t992 + t764 * t797;
t701 = t725 * t771 + t743 * t768;
t967 = t701 * t1025;
t966 = t701 * t1022;
t702 = -t725 * t768 + t743 * t771;
t965 = t702 * t1025;
t964 = t702 * t1022;
t726 = t747 * t800 - t792 * t806;
t744 = pkin(2) * t991 + t764 * t799;
t703 = t726 * t772 + t744 * t769;
t963 = t703 * t1024;
t962 = t703 * t1021;
t704 = -t726 * t769 + t744 * t772;
t961 = t704 * t1024;
t960 = t704 * t1021;
t727 = t748 * t802 - t793 * t808;
t745 = pkin(2) * t990 + t764 * t801;
t705 = t727 * t773 + t745 * t770;
t959 = t705 * t1023;
t958 = t705 * t1020;
t706 = -t727 * t770 + t745 * t773;
t957 = t706 * t1023;
t956 = t706 * t1020;
t955 = t713 * t1031;
t954 = t714 * t1029;
t953 = t715 * t1027;
t952 = t716 * t1031;
t951 = t717 * t1029;
t950 = t718 * t1027;
t949 = t734 * t1025;
t948 = t734 * t1022;
t947 = t735 * t1024;
t946 = t735 * t1021;
t945 = t736 * t1023;
t944 = t736 * t1020;
t943 = t755 * t1031;
t942 = t797 * t1031;
t941 = t737 * t1006;
t780 = t797 ^ 2;
t940 = t780 * t1030;
t939 = t757 * t1029;
t938 = t799 * t1029;
t937 = t739 * t1001;
t781 = t799 ^ 2;
t936 = t781 * t1028;
t935 = t759 * t1027;
t934 = t801 * t1027;
t933 = t741 * t996;
t782 = t801 ^ 2;
t932 = t782 * t1026;
t931 = t749 * t1006;
t930 = t750 * t1001;
t929 = t751 * t996;
t928 = t752 * t1006;
t927 = t753 * t1001;
t926 = t754 * t996;
t925 = t755 * t1008;
t924 = t757 * t1003;
t923 = t759 * t998;
t922 = t775 * t986;
t921 = t777 * t985;
t920 = t779 * t984;
t919 = t755 * t982;
t918 = t755 * t980;
t917 = t757 * t982;
t916 = t757 * t980;
t915 = t759 * t982;
t914 = t759 * t980;
t913 = t783 * t982;
t912 = t785 * t982;
t911 = t787 * t982;
t910 = t1027 * t1058;
t909 = t1029 * t1057;
t908 = t1031 * t1056;
t904 = t768 * t976;
t903 = t771 * t976;
t902 = t769 * t975;
t901 = t772 * t975;
t900 = t770 * t974;
t899 = t773 * t974;
t898 = qJ(3,1) * t934;
t897 = qJ(3,2) * t938;
t896 = qJ(3,3) * t942;
t895 = t1030 * t1055;
t894 = t1028 * t1054;
t893 = t1026 * t1053;
t892 = t737 * t949;
t891 = t737 * t948;
t890 = t739 * t947;
t889 = t739 * t946;
t888 = t741 * t945;
t887 = t741 * t944;
t886 = t755 * t942;
t885 = t803 * t943;
t884 = t780 * t941;
t883 = t738 * t922;
t882 = t757 * t938;
t881 = t805 * t939;
t880 = t781 * t937;
t879 = t740 * t921;
t878 = t759 * t934;
t877 = t807 * t935;
t876 = t782 * t933;
t875 = t742 * t920;
t874 = t797 * t925;
t873 = t803 * t925;
t872 = t799 * t924;
t871 = t805 * t924;
t870 = t801 * t923;
t869 = t807 * t923;
t868 = t942 * t1052;
t867 = t938 * t1052;
t866 = t934 * t1052;
t865 = t933 * t1058;
t864 = t937 * t1057;
t863 = t941 * t1056;
t862 = t768 * t886;
t861 = t771 * t886;
t860 = t737 * t804 * t922;
t859 = t769 * t882;
t858 = t772 * t882;
t857 = t739 * t806 * t921;
t856 = t770 * t878;
t855 = t773 * t878;
t854 = t741 * t808 * t920;
t853 = t768 * t874;
t852 = t771 * t874;
t851 = t769 * t872;
t850 = t772 * t872;
t849 = t770 * t870;
t848 = t773 * t870;
t847 = t907 * t737;
t846 = t906 * t739;
t845 = t905 * t741;
t844 = t913 * t1031;
t843 = t912 * t1029;
t842 = t911 * t1027;
t841 = qJ(3,1) * t754;
t840 = qJ(3,2) * t753;
t839 = qJ(3,3) * t752;
t838 = t759 * t841;
t837 = t778 * t841;
t836 = qJ(3,1) * t778 * t751;
t835 = t757 * t840;
t834 = t776 * t840;
t833 = qJ(3,2) * t776 * t750;
t832 = t755 * t839;
t831 = t774 * t839;
t830 = qJ(3,3) * t774 * t749;
t829 = t741 * t837;
t828 = t741 * t836;
t827 = t739 * t834;
t826 = t739 * t833;
t825 = t737 * t831;
t824 = t737 * t830;
t823 = (t713 * t768 + t716 * t771) * t943;
t822 = (t714 * t769 + t717 * t772) * t939;
t821 = (t715 * t770 + t718 * t773) * t935;
t820 = (t783 * t790 - t789 * t986) * t983;
t819 = (t785 * t790 - t789 * t985) * t983;
t818 = (t787 * t790 - t789 * t984) * t983;
t817 = t737 * t820;
t816 = t739 * t819;
t815 = t741 * t818;
t696 = t849 + t851 + t853;
t698 = t848 + t850 + t852;
t814 = t713 * t861 + t714 * t858 + t715 * t855;
t813 = t716 * t862 + t717 * t859 + t718 * t856;
t760 = 0.1e1 / t763 ^ 2;
t758 = 0.1e1 / t762 ^ 2;
t756 = 0.1e1 / t761 ^ 2;
t733 = t808 * t837;
t732 = t808 * t836;
t731 = t806 * t834;
t730 = t806 * t833;
t729 = t804 * t831;
t728 = t804 * t830;
t724 = (-t808 * t1048 + t736) * t778;
t723 = (-t806 * t1049 + t735) * t776;
t722 = (-t804 * t1050 + t734) * t774;
t712 = t718 ^ 2;
t711 = t717 ^ 2;
t710 = t716 ^ 2;
t709 = t715 ^ 2;
t708 = t714 ^ 2;
t707 = t713 ^ 2;
t700 = t756 * t768 * t771 + t758 * t769 * t772 + t760 * t770 * t773;
t699 = t771 * t873 + t772 * t871 + t773 * t869;
t697 = t768 * t873 + t769 * t871 + t770 * t869;
t695 = (t808 * t911 + (t801 * t808 * t981 + t736 * t789) * t807 + t736 * t987) * t778;
t694 = (t806 * t912 + (t799 * t806 * t981 + t735 * t789) * t805 + t735 * t988) * t776;
t693 = (t804 * t913 + (t797 * t804 * t981 + t734 * t789) * t803 + t734 * t989) * t774;
t692 = (t751 * t736 + t808 * t818) * t778;
t691 = (t750 * t735 + t806 * t819) * t776;
t690 = (t749 * t734 + t804 * t820) * t774;
t689 = t715 * t910 - t899;
t688 = t714 * t909 - t901;
t687 = t713 * t908 - t903;
t686 = t718 * t910 - t900;
t685 = t717 * t909 - t902;
t684 = t716 * t908 - t904;
t683 = pkin(1) * t1015 - t718 * t898;
t682 = pkin(1) * t1014 - t715 * t898;
t681 = pkin(1) * t1017 - t717 * t897;
t680 = pkin(1) * t1016 - t714 * t897;
t679 = pkin(1) * t1019 - t716 * t896;
t678 = pkin(1) * t1018 - t713 * t896;
t677 = (-t718 * t977 + t705) * t778;
t676 = (-t715 * t977 + t706) * t778;
t675 = (-t717 * t978 + t703) * t776;
t674 = (-t714 * t978 + t704) * t776;
t673 = (-t716 * t979 + t701) * t774;
t672 = (-t713 * t979 + t702) * t774;
t671 = -t718 * t829 + t770 * t914;
t670 = -t715 * t829 + t773 * t914;
t669 = t718 * t828 + t770 * t915;
t668 = t715 * t828 + t773 * t915;
t667 = -t717 * t827 + t769 * t916;
t666 = -t714 * t827 + t772 * t916;
t665 = t717 * t826 + t769 * t917;
t664 = t714 * t826 + t772 * t917;
t663 = -t716 * t825 + t768 * t918;
t662 = -t713 * t825 + t771 * t918;
t661 = t716 * t824 + t768 * t919;
t660 = t713 * t824 + t771 * t919;
t659 = t716 * t941 + t717 * t937 + t718 * t933;
t658 = t713 * t941 + t714 * t937 + t715 * t933;
t657 = -qJ(3,1) * t900 + (-t705 * t1048 + t718 * t845) * t778;
t656 = -qJ(3,1) * t899 + (-t706 * t1048 + t715 * t845) * t778;
t655 = -qJ(3,2) * t902 + (-t703 * t1049 + t717 * t846) * t776;
t654 = -qJ(3,2) * t901 + (-t704 * t1049 + t714 * t846) * t776;
t653 = -qJ(3,3) * t904 + (-t701 * t1050 + t716 * t847) * t774;
t652 = -qJ(3,3) * t903 + (-t702 * t1050 + t713 * t847) * t774;
t651 = t716 * t884 + t717 * t880 + t718 * t876;
t650 = t713 * t884 + t714 * t880 + t715 * t876;
t649 = 0.2e1 * t716 * t860 + 0.2e1 * t717 * t857 + 0.2e1 * t718 * t854;
t648 = 0.2e1 * t713 * t860 + 0.2e1 * t714 * t857 + 0.2e1 * t715 * t854;
t647 = t718 * t842 + ((-qJ(3,1) * t1015 + t718 * t866) * t790 + t705 * t1000) * t807 + t801 * (t705 * t999 + t770 * t970);
t646 = t715 * t842 + ((-qJ(3,1) * t1014 + t715 * t866) * t790 + t706 * t1000) * t807 + t801 * (t706 * t999 + t773 * t970);
t645 = t717 * t843 + ((-qJ(3,2) * t1017 + t717 * t867) * t790 + t703 * t1005) * t805 + t799 * (t703 * t1004 + t769 * t969);
t644 = t714 * t843 + ((-qJ(3,2) * t1016 + t714 * t867) * t790 + t704 * t1005) * t805 + t799 * (t704 * t1004 + t772 * t969);
t643 = t716 * t844 + ((-qJ(3,3) * t1019 + t716 * t868) * t790 + t701 * t1010) * t803 + t797 * (t701 * t1009 + t768 * t968);
t642 = t713 * t844 + ((-qJ(3,3) * t1018 + t713 * t868) * t790 + t702 * t1010) * t803 + t797 * (t702 * t1009 + t771 * t968);
t641 = -t770 * t838 + (t751 * t705 + t718 * t815) * t778;
t640 = -t773 * t838 + (t751 * t706 + t715 * t815) * t778;
t639 = -t769 * t835 + (t750 * t703 + t717 * t816) * t776;
t638 = -t772 * t835 + (t750 * t704 + t714 * t816) * t776;
t637 = -t768 * t832 + (t749 * t701 + t716 * t817) * t774;
t636 = -t771 * t832 + (t749 * t702 + t713 * t817) * t774;
t635 = t893 + t894 + t895;
t634 = t780 * t895 + t781 * t894 + t782 * t893;
t633 = 0.2e1 * t875 * t1053 + 0.2e1 * t879 * t1054 + 0.2e1 * t883 * t1055;
t632 = t803 * t823 + t805 * t822 + t807 * t821;
t631 = t797 * t823 + t799 * t822 + t801 * t821;
t1 = [t712 * t1026 + t711 * t1028 + t710 * t1030, 0, 0, t710 * t940 + t711 * t936 + t712 * t932, 0.2e1 * t710 * t883 + 0.2e1 * t711 * t879 + 0.2e1 * t712 * t875, 0.2e1 * t813, 0.2e1 * t716 * t768 * t885 + 0.2e1 * t717 * t769 * t881 + 0.2e1 * t718 * t770 * t877, t756 * t768 ^ 2 + t758 * t769 ^ 2 + t760 * t770 ^ 2, 0, 0, t663 * t1019 + t667 * t1017 + t671 * t1015 + (t641 * t778 + t959) * t1036 + (t639 * t776 + t963) * t1038 + (t637 * t774 + t967) * t1040, t661 * t1019 + t665 * t1017 + t669 * t1015 + (t647 * t778 + t958) * t1036 + (t645 * t776 + t962) * t1038 + (t643 * t774 + t966) * t1040, -t813 * pkin(1) + t684 * t952 + t685 * t951 + t686 * t950, (t657 * t1036 + t677 * t705) * t778 + (t655 * t1038 + t675 * t703) * t776 + (t653 * t1040 + t673 * t701) * t774 + (t683 * t1015 + t681 * t1017 + t679 * t1019) * pkin(1), 1; t635, 0, 0, t634, t633, t631, t632, t700, 0, 0, t662 * t1019 + t666 * t1017 + t670 * t1015 + (t640 * t1035 + t715 * t959) * t741 + (t638 * t1037 + t714 * t963) * t739 + (t636 * t1039 + t713 * t967) * t737, t660 * t1019 + t664 * t1017 + t668 * t1015 + (t646 * t1035 + t715 * t958) * t741 + (t644 * t1037 + t714 * t962) * t739 + (t642 * t1039 + t713 * t966) * t737, t687 * t952 + t688 * t951 + t689 * t950 + (-t713 * t862 - t714 * t859 - t715 * t856) * pkin(1), (t656 * t1036 + t676 * t705) * t778 + (t654 * t1038 + t674 * t703) * t776 + (t652 * t1040 + t672 * t701) * t774 + (t682 * t1015 + t680 * t1017 + t678 * t1019) * pkin(1), 0; t659, 0, 0, t651, t649, t696, t697, 0, 0, 0, -t733 * t1015 - t731 * t1017 - t729 * t1019 + t690 * t952 + t691 * t951 + t692 * t950 + t701 * t931 + t703 * t930 + t705 * t929, t732 * t1015 + t730 * t1017 + t728 * t1019 + t693 * t952 + t694 * t951 + t695 * t950 + t701 * t928 + t703 * t927 + t705 * t926, -t696 * pkin(1) + t716 * t863 + t717 * t864 + t718 * t865, (t718 * t1032 + t705 * t724) * t778 + (t717 * t1033 + t703 * t723) * t776 + (t716 * t1034 + t701 * t722) * t774 + (-qJ(3,1) * t849 - qJ(3,2) * t851 - qJ(3,3) * t853) * pkin(1), 0; t635, 0, 0, t634, t633, t631, t632, t700, 0, 0, t663 * t1018 + t667 * t1016 + t671 * t1014 + (t641 * t1041 + t718 * t957) * t741 + (t639 * t1043 + t717 * t961) * t739 + (t637 * t1045 + t716 * t965) * t737, t661 * t1018 + t665 * t1016 + t669 * t1014 + (t647 * t1041 + t718 * t956) * t741 + (t645 * t1043 + t717 * t960) * t739 + (t643 * t1045 + t716 * t964) * t737, t684 * t955 + t685 * t954 + t686 * t953 + (-t716 * t861 - t717 * t858 - t718 * t855) * pkin(1), (t657 * t1042 + t677 * t706) * t778 + (t655 * t1044 + t675 * t704) * t776 + (t653 * t1046 + t673 * t702) * t774 + (t683 * t1014 + t681 * t1016 + t679 * t1018) * pkin(1), 0; t709 * t1026 + t708 * t1028 + t707 * t1030, 0, 0, t707 * t940 + t708 * t936 + t709 * t932, 0.2e1 * t707 * t883 + 0.2e1 * t708 * t879 + 0.2e1 * t709 * t875, 0.2e1 * t814, 0.2e1 * t713 * t771 * t885 + 0.2e1 * t714 * t772 * t881 + 0.2e1 * t715 * t773 * t877, t756 * t771 ^ 2 + t758 * t772 ^ 2 + t760 * t773 ^ 2, 0, 0, t662 * t1018 + t666 * t1016 + t670 * t1014 + (t640 * t778 + t957) * t1042 + (t638 * t776 + t961) * t1044 + (t636 * t774 + t965) * t1046, t660 * t1018 + t664 * t1016 + t668 * t1014 + (t646 * t778 + t956) * t1042 + (t644 * t776 + t960) * t1044 + (t642 * t774 + t964) * t1046, -t814 * pkin(1) + t687 * t955 + t688 * t954 + t689 * t953, (t656 * t1042 + t676 * t706) * t778 + (t654 * t1044 + t674 * t704) * t776 + (t652 * t1046 + t672 * t702) * t774 + (t682 * t1014 + t680 * t1016 + t678 * t1018) * pkin(1), 1; t658, 0, 0, t650, t648, t698, t699, 0, 0, 0, -t733 * t1014 - t731 * t1016 - t729 * t1018 + t690 * t955 + t691 * t954 + t692 * t953 + t702 * t931 + t704 * t930 + t706 * t929, t732 * t1014 + t730 * t1016 + t728 * t1018 + t693 * t955 + t694 * t954 + t695 * t953 + t702 * t928 + t704 * t927 + t706 * t926, -t698 * pkin(1) + t713 * t863 + t714 * t864 + t715 * t865, (t715 * t1032 + t706 * t724) * t778 + (t714 * t1033 + t704 * t723) * t776 + (t713 * t1034 + t702 * t722) * t774 + (-qJ(3,1) * t848 - qJ(3,2) * t850 - qJ(3,3) * t852) * pkin(1), 0; t659, 0, 0, t651, t649, t696, t697, 0, 0, 0, t639 * t1003 + t637 * t1008 + t641 * t998 + t716 * t892 + t717 * t890 + t718 * t888, t645 * t1003 + t643 * t1008 + t647 * t998 + t716 * t891 + t717 * t889 + t718 * t887, t685 * t1003 + t684 * t1008 + t686 * t998, (t657 * t808 + t677 * t736) * t778 + (t655 * t806 + t675 * t735) * t776 + (t653 * t804 + t673 * t734) * t774, 0; t658, 0, 0, t650, t648, t698, t699, 0, 0, 0, t638 * t1003 + t636 * t1008 + t640 * t998 + t713 * t892 + t714 * t890 + t715 * t888, t644 * t1003 + t642 * t1008 + t646 * t998 + t713 * t891 + t714 * t889 + t715 * t887, t688 * t1003 + t687 * t1008 + t689 * t998, (t656 * t808 + t676 * t736) * t778 + (t654 * t806 + t674 * t735) * t776 + (t652 * t804 + t672 * t734) * t774, 0; t1002 + t997 + t1007, 0, 0, t781 * t1002 + t780 * t1007 + t782 * t997, 0.2e1 * t784 * t922 + 0.2e1 * t786 * t921 + 0.2e1 * t788 * t920, 0, 0, 0, 0, 0, (t692 * t778 + t945) * t808 + (t691 * t776 + t947) * t806 + (t690 * t774 + t949) * t804, (t695 * t778 + t944) * t808 + (t694 * t776 + t946) * t806 + (t693 * t774 + t948) * t804, 0.2e1 * qJ(3,1) * t997 + 0.2e1 * qJ(3,2) * t1002 + 0.2e1 * qJ(3,3) * t1007, (t721 * t808 + t724 * t736) * t778 + (t720 * t806 + t723 * t735) * t776 + (t719 * t804 + t722 * t734) * t774, 1;];
tau_reg  = t1;
