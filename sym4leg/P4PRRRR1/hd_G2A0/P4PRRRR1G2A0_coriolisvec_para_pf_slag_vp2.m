% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:44
% EndTime: 2020-08-07 10:58:50
% DurationCPUTime: 6.41s
% Computational Cost: add. (5736->360), mult. (11772->758), div. (5748->13), fcn. (12420->26), ass. (0->323)
t832 = cos(qJ(3,4));
t813 = t832 ^ 2;
t1034 = 0.2e1 * t813;
t856 = xP(4);
t810 = sin(t856);
t811 = cos(t856);
t859 = koppelP(4,2);
t863 = koppelP(4,1);
t782 = t810 * t863 + t811 * t859;
t786 = -t810 * t859 + t811 * t863;
t836 = legFrame(4,2);
t802 = sin(t836);
t806 = cos(t836);
t852 = xDP(4);
t854 = xDP(2);
t855 = xDP(1);
t899 = (t786 * t852 + t854) * t802 - (-t782 * t852 + t855) * t806;
t814 = 0.1e1 / t832;
t867 = 0.1e1 / pkin(2);
t977 = t814 * t867;
t734 = t899 * t977;
t830 = sin(qJ(3,4));
t1007 = t734 * t830;
t831 = sin(qJ(2,4));
t812 = 0.1e1 / t831;
t853 = xDP(3);
t815 = 0.1e1 / t832 ^ 2;
t833 = cos(qJ(2,4));
t937 = t815 * t830 * t833;
t707 = (-t814 * t853 - t899 * t937) * t867 * t812;
t1021 = t707 * t833;
t966 = t831 * t832;
t662 = ((-t831 * t1007 + t832 * t1021) * t814 * t707 + (-t707 * t830 * t966 + t734 * t833) * t815 * t734) * t812;
t730 = t734 ^ 2;
t1011 = t730 * t814;
t706 = t707 ^ 2;
t682 = (-t706 * t832 - t1011) * t812 * pkin(2);
t790 = Ifges(3,5) * t830 + Ifges(3,6) * t832;
t955 = t830 * t1011;
t834 = Ifges(3,1) - Ifges(3,2);
t967 = t830 * t834;
t791 = mrSges(3,1) * t830 + mrSges(3,2) * t832;
t994 = t791 * t831;
t1051 = Ifges(3,3) * t955 - t662 * t790 + t682 * t994 - t706 * (Ifges(3,4) * t1034 + t832 * t967 - Ifges(3,4));
t862 = koppelP(1,2);
t866 = koppelP(1,1);
t785 = t810 * t866 + t811 * t862;
t789 = -t810 * t862 + t811 * t866;
t839 = legFrame(1,2);
t805 = sin(t839);
t809 = cos(t839);
t896 = (t789 * t852 + t854) * t805 - (-t785 * t852 + t855) * t809;
t850 = cos(qJ(3,1));
t827 = 0.1e1 / t850;
t968 = t827 * t867;
t737 = t896 * t968;
t844 = sin(qJ(3,1));
t1004 = t737 * t844;
t845 = sin(qJ(2,1));
t819 = 0.1e1 / t845;
t828 = 0.1e1 / t850 ^ 2;
t851 = cos(qJ(2,1));
t928 = t828 * t844 * t851;
t713 = (-t827 * t853 - t896 * t928) * t867 * t819;
t1013 = t713 * t844;
t1029 = -Ifges(3,6) / 0.2e1;
t1030 = Ifges(3,5) / 0.2e1;
t826 = t850 ^ 2;
t1031 = 0.2e1 * t826;
t1012 = t713 * t851;
t960 = t845 * t850;
t665 = ((-t845 * t1004 + t850 * t1012) * t827 * t713 + (-t960 * t1013 + t737 * t851) * t828 * t737) * t819;
t733 = t737 ^ 2;
t1008 = t733 * t827;
t710 = t713 ^ 2;
t685 = (-t710 * t850 - t1008) * t819 * pkin(2);
t835 = mrSges(2,2) - mrSges(3,3);
t912 = -mrSges(3,1) * t850 + mrSges(3,2) * t844;
t757 = (mrSges(2,1) - t912) * t851 - t845 * t835;
t1027 = Ifges(3,1) + Ifges(2,3);
t1035 = 0.2e1 * Ifges(3,4);
t961 = t844 * t850;
t777 = t961 * t1035 - t826 * t834 + t1027;
t794 = Ifges(3,5) * t844 + Ifges(3,6) * t850;
t952 = t844 * t1008;
t1043 = -t665 * t777 - t685 * t757 + t794 * t952 + 0.2e1 * ((t834 * t1013 + t737 * t1030) * t850 + t1004 * t1029 + (t1031 - 0.1e1) * t713 * Ifges(3,4)) * t737;
t861 = koppelP(2,2);
t865 = koppelP(2,1);
t784 = t810 * t865 + t811 * t861;
t788 = -t810 * t861 + t811 * t865;
t838 = legFrame(2,2);
t804 = sin(t838);
t808 = cos(t838);
t897 = (t788 * t852 + t854) * t804 - (-t784 * t852 + t855) * t808;
t848 = cos(qJ(3,2));
t824 = 0.1e1 / t848;
t969 = t824 * t867;
t736 = t897 * t969;
t842 = sin(qJ(3,2));
t1005 = t736 * t842;
t843 = sin(qJ(2,2));
t818 = 0.1e1 / t843;
t825 = 0.1e1 / t848 ^ 2;
t849 = cos(qJ(2,2));
t929 = t825 * t842 * t849;
t712 = (-t824 * t853 - t897 * t929) * t867 * t818;
t1015 = t712 * t842;
t823 = t848 ^ 2;
t1032 = 0.2e1 * t823;
t1014 = t712 * t849;
t962 = t843 * t848;
t664 = ((-t843 * t1005 + t848 * t1014) * t824 * t712 + (-t962 * t1015 + t736 * t849) * t825 * t736) * t818;
t732 = t736 ^ 2;
t1009 = t732 * t824;
t709 = t712 ^ 2;
t684 = (-t709 * t848 - t1009) * t818 * pkin(2);
t913 = -mrSges(3,1) * t848 + mrSges(3,2) * t842;
t756 = (mrSges(2,1) - t913) * t849 - t843 * t835;
t963 = t842 * t848;
t776 = t963 * t1035 - t823 * t834 + t1027;
t793 = Ifges(3,5) * t842 + Ifges(3,6) * t848;
t953 = t842 * t1009;
t1042 = -t664 * t776 - t684 * t756 + t793 * t953 + 0.2e1 * ((t834 * t1015 + t736 * t1030) * t848 + t1005 * t1029 + (t1032 - 0.1e1) * t712 * Ifges(3,4)) * t736;
t860 = koppelP(3,2);
t864 = koppelP(3,1);
t783 = t810 * t864 + t811 * t860;
t787 = -t810 * t860 + t811 * t864;
t837 = legFrame(3,2);
t803 = sin(t837);
t807 = cos(t837);
t898 = (t787 * t852 + t854) * t803 - (-t783 * t852 + t855) * t807;
t846 = cos(qJ(3,3));
t821 = 0.1e1 / t846;
t970 = t821 * t867;
t735 = t898 * t970;
t840 = sin(qJ(3,3));
t1006 = t735 * t840;
t841 = sin(qJ(2,3));
t817 = 0.1e1 / t841;
t822 = 0.1e1 / t846 ^ 2;
t847 = cos(qJ(2,3));
t930 = t822 * t840 * t847;
t711 = (-t821 * t853 - t898 * t930) * t867 * t817;
t1017 = t711 * t840;
t820 = t846 ^ 2;
t1033 = 0.2e1 * t820;
t1016 = t711 * t847;
t964 = t841 * t846;
t663 = ((-t841 * t1006 + t846 * t1016) * t821 * t711 + (-t964 * t1017 + t735 * t847) * t822 * t735) * t817;
t731 = t735 ^ 2;
t1010 = t731 * t821;
t708 = t711 ^ 2;
t683 = (-t708 * t846 - t1010) * t817 * pkin(2);
t914 = -mrSges(3,1) * t846 + mrSges(3,2) * t840;
t755 = (mrSges(2,1) - t914) * t847 - t841 * t835;
t965 = t840 * t846;
t775 = t965 * t1035 - t820 * t834 + t1027;
t792 = Ifges(3,5) * t840 + Ifges(3,6) * t846;
t954 = t840 * t1010;
t1041 = -t663 * t775 - t683 * t755 + t792 * t954 + 0.2e1 * ((t834 * t1017 + t735 * t1030) * t846 + t1006 * t1029 + (t1033 - 0.1e1) * t711 * Ifges(3,4)) * t735;
t915 = -mrSges(3,1) * t832 + mrSges(3,2) * t830;
t754 = (mrSges(2,1) - t915) * t833 - t831 * t835;
t774 = t830 * t832 * t1035 - t813 * t834 + t1027;
t1040 = -t662 * t774 - t682 * t754 + t790 * t955 + 0.2e1 * ((t734 * t1030 + t707 * t967) * t832 + t1007 * t1029 + (t1034 - 0.1e1) * t707 * Ifges(3,4)) * t734;
t795 = mrSges(3,1) * t840 + mrSges(3,2) * t846;
t992 = t795 * t841;
t1050 = Ifges(3,3) * t954 - t663 * t792 + t683 * t992 - t708 * (Ifges(3,4) * t1033 + t834 * t965 - Ifges(3,4));
t796 = mrSges(3,1) * t842 + mrSges(3,2) * t848;
t990 = t796 * t843;
t1049 = Ifges(3,3) * t953 - t664 * t793 + t684 * t990 - t709 * (Ifges(3,4) * t1032 + t834 * t963 - Ifges(3,4));
t797 = mrSges(3,1) * t844 + mrSges(3,2) * t850;
t988 = t797 * t845;
t1048 = Ifges(3,3) * t952 - t665 * t794 + t685 * t988 - t710 * (Ifges(3,4) * t1031 + t834 * t961 - Ifges(3,4));
t1028 = t835 / 0.2e1;
t816 = m(1) + m(2) + m(3);
t943 = t814 * t994;
t1047 = -t730 * t830 * t943 - t662 * t754 - t682 * t816 + (-mrSges(2,1) * t706 + t915 * (t706 + t730)) * t831 - 0.2e1 * (t707 * t1028 + t791 * t734) * t1021;
t942 = t821 * t992;
t1046 = -t731 * t840 * t942 - t663 * t755 - t683 * t816 + (-mrSges(2,1) * t708 + t914 * (t708 + t731)) * t841 - 0.2e1 * (t711 * t1028 + t795 * t735) * t1016;
t941 = t824 * t990;
t1045 = -t732 * t842 * t941 - t664 * t756 - t684 * t816 + (-mrSges(2,1) * t709 + t913 * (t709 + t732)) * t843 - 0.2e1 * (t712 * t1028 + t796 * t736) * t1014;
t940 = t827 * t988;
t1044 = -t733 * t844 * t940 - t665 * t757 - t685 * t816 + (-mrSges(2,1) * t710 + t912 * (t710 + t733)) * t845 - 0.2e1 * (t713 * t1028 + t797 * t737) * t1012;
t916 = t867 * t928;
t900 = t819 * t916;
t1039 = -t1043 * t900 + t1048 * t968;
t917 = t867 * t929;
t901 = t818 * t917;
t1038 = -t1042 * t901 + t1049 * t969;
t918 = t867 * t930;
t902 = t817 * t918;
t1037 = -t1041 * t902 + t1050 * t970;
t922 = t867 * t937;
t903 = t812 * t922;
t1036 = -t1040 * t903 + t1051 * t977;
t829 = t852 ^ 2;
t766 = t802 * t966 - t806 * t830;
t1003 = t766 * t814;
t767 = t802 * t830 + t806 * t966;
t1002 = t767 * t814;
t768 = t803 * t964 - t807 * t840;
t1001 = t768 * t821;
t769 = t804 * t962 - t808 * t842;
t1000 = t769 * t824;
t770 = t805 * t960 - t809 * t844;
t999 = t770 * t827;
t771 = t803 * t840 + t807 * t964;
t998 = t771 * t821;
t772 = t804 * t842 + t808 * t962;
t997 = t772 * t824;
t773 = t805 * t844 + t809 * t960;
t996 = t773 * t827;
t995 = t791 * t814;
t993 = t795 * t821;
t991 = t796 * t824;
t989 = t797 * t827;
t987 = t802 * t867;
t986 = t803 * t867;
t985 = t804 * t867;
t984 = t805 * t867;
t983 = t806 * t867;
t982 = t807 * t867;
t981 = t808 * t867;
t980 = t809 * t867;
t979 = t812 * t814;
t976 = t817 * t821;
t974 = t818 * t824;
t972 = t819 * t827;
t951 = t766 * t979;
t950 = t767 * t979;
t949 = t768 * t976;
t948 = t769 * t974;
t947 = t770 * t972;
t946 = t771 * t976;
t945 = t772 * t974;
t944 = t773 * t972;
t939 = t816 * t979;
t938 = t812 * t977;
t936 = t816 * t976;
t935 = t816 * t974;
t934 = t816 * t972;
t933 = t817 * t970;
t932 = t818 * t969;
t931 = t819 * t968;
t923 = t812 * t937;
t921 = t817 * t930;
t920 = t818 * t929;
t919 = t819 * t928;
t911 = t802 * t922;
t910 = t803 * t918;
t909 = t804 * t917;
t908 = t805 * t916;
t907 = t806 * t922;
t906 = t807 * t918;
t905 = t808 * t917;
t904 = t809 * t916;
t895 = t782 * t806 + t786 * t802;
t894 = t783 * t807 + t787 * t803;
t893 = t784 * t808 + t788 * t804;
t892 = t785 * t809 + t789 * t805;
t883 = (-t782 * t802 + t786 * t806) * t977;
t882 = (-t783 * t803 + t787 * t807) * t970;
t881 = (-t784 * t804 + t788 * t808) * t969;
t880 = (-t785 * t805 + t789 * t809) * t968;
t879 = -Ifges(3,3) * t814 + t790 * t923;
t878 = -Ifges(3,3) * t821 + t792 * t921;
t877 = -Ifges(3,3) * t824 + t793 * t920;
t876 = -Ifges(3,3) * t827 + t794 * t919;
t875 = t774 * t923 - t790 * t814;
t874 = t775 * t921 - t792 * t821;
t873 = t776 * t920 - t793 * t824;
t872 = t777 * t919 - t794 * t827;
t871 = t754 * t923 + t943;
t870 = t755 * t921 + t942;
t869 = t756 * t920 + t941;
t868 = t757 * t919 + t940;
t858 = mrSges(4,1);
t857 = mrSges(4,2);
t749 = (-t757 * t968 + t816 * t851) * t819;
t748 = (-t756 * t969 + t816 * t849) * t818;
t747 = (-t755 * t970 + t816 * t847) * t817;
t746 = (-t754 * t977 + t816 * t833) * t812;
t745 = t892 * t968;
t744 = t893 * t969;
t743 = t894 * t970;
t742 = t895 * t977;
t741 = (t757 * t851 - t777 * t968) * t819;
t740 = (t756 * t849 - t776 * t969) * t818;
t739 = (t755 * t847 - t775 * t970) * t817;
t738 = (t754 * t833 - t774 * t977) * t812;
t729 = t892 * t900;
t728 = t893 * t901;
t727 = t894 * t902;
t726 = t895 * t903;
t717 = (-t770 * t785 + t773 * t789) * t972;
t716 = (-t769 * t784 + t772 * t788) * t974;
t715 = (-t768 * t783 + t771 * t787) * t976;
t714 = (-t766 * t782 + t767 * t786) * t979;
t705 = t773 * t934 - t868 * t984;
t704 = t772 * t935 - t869 * t985;
t703 = t771 * t936 - t870 * t986;
t702 = t770 * t934 + t868 * t980;
t701 = t769 * t935 + t869 * t981;
t700 = t768 * t936 + t870 * t982;
t699 = t767 * t939 - t871 * t987;
t698 = t766 * t939 + t871 * t983;
t697 = t757 * t944 - t872 * t984;
t696 = t756 * t945 - t873 * t985;
t695 = t755 * t946 - t874 * t986;
t694 = t757 * t947 + t872 * t980;
t693 = t756 * t948 + t873 * t981;
t692 = t755 * t949 + t874 * t982;
t691 = t754 * t950 - t875 * t987;
t690 = t754 * t951 + t875 * t983;
t677 = t717 * t816 - t729 * t757 - t745 * t988;
t676 = t716 * t816 - t728 * t756 - t744 * t990;
t675 = t715 * t816 - t727 * t755 - t743 * t992;
t674 = t714 * t816 - t726 * t754 - t742 * t994;
t673 = t717 * t757 - t729 * t777 + t745 * t794;
t672 = t716 * t756 - t728 * t776 + t744 * t793;
t671 = t715 * t755 - t727 * t775 + t743 * t792;
t670 = t714 * t754 - t726 * t774 + t742 * t790;
t1 = [t1047 * t951 + t1046 * t949 + t1045 * t948 + t1044 * t947 - t1039 * t809 - t1038 * t808 - t1037 * t807 - t1036 * t806 + (t810 * t857 - t811 * t858 + (-t766 * t995 + t879 * t983) * t883 + (-(t698 * t1003 + t690 * t907) * t786 - (t698 * t1002 - t690 * t911) * t782) * t812 + (-t769 * t991 + t877 * t981) * t881 + (-(t701 * t1000 + t693 * t905) * t788 - (-t693 * t909 + t701 * t997) * t784) * t818 + (-t768 * t993 + t878 * t982) * t882 + (-(t700 * t1001 + t692 * t906) * t787 - (-t692 * t910 + t700 * t998) * t783) * t817 + (-t770 * t989 + t876 * t980) * t880 + (-(t694 * t904 + t702 * t999) * t789 - (-t694 * t908 + t702 * t996) * t785) * t819) * t829; t1047 * t950 + t1046 * t946 + t1045 * t945 + t1044 * t944 + t1039 * t805 + t1038 * t804 + t1037 * t803 + t1036 * t802 + (-t810 * t858 - t811 * t857 + (-t767 * t995 - t879 * t987) * t883 + (-(t699 * t1003 + t691 * t907) * t786 - (t699 * t1002 - t691 * t911) * t782) * t812 + (-t772 * t991 - t877 * t985) * t881 + (-(t704 * t1000 + t696 * t905) * t788 - (-t696 * t909 + t704 * t997) * t784) * t818 + (-t771 * t993 - t878 * t986) * t882 + (-(t703 * t1001 + t695 * t906) * t787 - (-t695 * t910 + t703 * t998) * t783) * t817 + (-t773 * t989 - t876 * t984) * t880 + (-(t697 * t904 + t705 * t999) * t789 - (-t697 * t908 + t705 * t996) * t785) * t819) * t829; t1047 * t812 * t833 + t1046 * t817 * t847 + t1045 * t818 * t849 + t1044 * t819 * t851 - t1040 * t938 - t1041 * t933 - t1042 * t932 - t1043 * t931 + ((-t794 * t931 - t851 * t797) * t880 + (-(t741 * t904 + t749 * t999) * t789 - (-t741 * t908 + t749 * t996) * t785) * t819 + (-t793 * t932 - t849 * t796) * t881 + (-(t748 * t1000 + t740 * t905) * t788 - (-t740 * t909 + t748 * t997) * t784) * t818 + (-t792 * t933 - t847 * t795) * t882 + (-(t747 * t1001 + t739 * t906) * t787 - (-t739 * t910 + t747 * t998) * t783) * t817 + (-t790 * t938 - t833 * t791) * t883 + (-(t746 * t1003 + t738 * t907) * t786 - (t746 * t1002 - t738 * t911) * t782) * t812) * t829; t1048 * t745 + t1049 * t744 + t1050 * t743 + t1051 * t742 - t1043 * t729 - t1042 * t728 - t1041 * t727 - t1040 * t726 + t1044 * t717 + t1045 * t716 + t1046 * t715 + t1047 * t714 + ((Ifges(3,3) * t742 - t714 * t994 - t726 * t790) * t883 + (-(t674 * t1003 + t670 * t907) * t786 - (t674 * t1002 - t670 * t911) * t782) * t812 + (Ifges(3,3) * t745 - t717 * t988 - t729 * t794) * t880 + (-(t673 * t904 + t677 * t999) * t789 - (-t673 * t908 + t677 * t996) * t785) * t819 + (Ifges(3,3) * t744 - t716 * t990 - t728 * t793) * t881 + (-(t676 * t1000 + t672 * t905) * t788 - (-t672 * t909 + t676 * t997) * t784) * t818 + (Ifges(3,3) * t743 - t715 * t992 - t727 * t792) * t882 + (-(t675 * t1001 + t671 * t906) * t787 - (-t671 * t910 + t675 * t998) * t783) * t817) * t829;];
taucX  = t1;
