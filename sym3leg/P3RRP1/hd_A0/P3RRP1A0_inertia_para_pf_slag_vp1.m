% Calculate inertia matrix for parallel robot
% P3RRP1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m [4x1]
%   mass of all robot links (including platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 18:08
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MX = P3RRP1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:59:31
% EndTime: 2018-12-20 17:59:35
% DurationCPUTime: 3.87s
% Computational Cost: add. (24448->471), mult. (51855->673), div. (432->3), fcn. (21590->20), ass. (0->297)
t900 = (qJ(3,3) ^ 2);
t912 = pkin(1) ^ 2;
t891 = cos(qJ(1,3));
t875 = t891 ^ 2;
t953 = t875 * t912;
t1010 = (2 * t900) + t912 - t953;
t901 = (qJ(3,2) ^ 2);
t892 = cos(qJ(1,2));
t876 = t892 ^ 2;
t952 = t876 * t912;
t1009 = (2 * t901) + t912 - t952;
t902 = (qJ(3,1) ^ 2);
t893 = cos(qJ(1,1));
t877 = t893 ^ 2;
t951 = t877 * t912;
t1008 = (2 * t902) + t912 - t951;
t1007 = 2 * pkin(1);
t1006 = 2 * rSges(3,3);
t885 = legFrame(3,3);
t854 = sin(t885);
t834 = t854 * qJ(3,3);
t831 = pkin(2) * t834;
t1005 = -0.2e1 * t831;
t886 = legFrame(2,3);
t855 = sin(t886);
t835 = t855 * qJ(3,2);
t832 = pkin(2) * t835;
t1004 = -0.2e1 * t832;
t887 = legFrame(1,3);
t856 = sin(t887);
t836 = t856 * qJ(3,1);
t833 = pkin(2) * t836;
t1003 = -0.2e1 * t833;
t888 = sin(qJ(1,3));
t1002 = -0.2e1 * t888;
t889 = sin(qJ(1,2));
t1001 = -0.2e1 * t889;
t890 = sin(qJ(1,1));
t1000 = -0.2e1 * t890;
t999 = -0.2e1 * t891;
t998 = -0.2e1 * t892;
t997 = -0.2e1 * t893;
t894 = pkin(2) + rSges(3,1);
t878 = qJ(1,3) + qJ(2,3);
t851 = cos(t878);
t842 = t851 ^ 2;
t911 = pkin(2) ^ 2;
t860 = -t900 + t911;
t845 = 0.1e1 + t860;
t848 = sin(t878);
t881 = 0.2e1 + t912;
t950 = t888 * t891;
t987 = pkin(2) * qJ(3,3);
t913 = -t845 * t875 + 0.2e1 * t950 * t987;
t947 = -0.1e1 / 0.2e1 - t911 / 0.2e1;
t965 = t848 * qJ(3,3);
t975 = qJ(3,3) * t888;
t747 = 0.1e1 / (((t845 * t950 + (0.2e1 * t875 - 0.1e1) * t987) * pkin(1) * t848 + t975) * t851 * t1007 + pkin(1) * t965 * t999 - t900 * t881 + (0.2e1 * (t900 / 0.2e1 - t913 + t947) * t842 + t913) * t912);
t996 = m(3) * t747;
t879 = qJ(1,2) + qJ(2,2);
t852 = cos(t879);
t843 = t852 ^ 2;
t861 = -t901 + t911;
t846 = 0.1e1 + t861;
t849 = sin(t879);
t949 = t889 * t892;
t988 = pkin(2) * qJ(3,2);
t914 = -t846 * t876 + 0.2e1 * t949 * t988;
t963 = t849 * qJ(3,2);
t977 = qJ(3,2) * t889;
t748 = 0.1e1 / (((t846 * t949 + (0.2e1 * t876 - 0.1e1) * t988) * pkin(1) * t849 + t977) * t852 * t1007 + pkin(1) * t963 * t998 - t901 * t881 + (0.2e1 * (t901 / 0.2e1 - t914 + t947) * t843 + t914) * t912);
t995 = m(3) * t748;
t880 = qJ(1,1) + qJ(2,1);
t853 = cos(t880);
t844 = t853 ^ 2;
t862 = -t902 + t911;
t847 = 0.1e1 + t862;
t850 = sin(t880);
t948 = t890 * t893;
t989 = pkin(2) * qJ(3,1);
t915 = -t847 * t877 + 0.2e1 * t948 * t989;
t961 = t850 * qJ(3,1);
t979 = qJ(3,1) * t890;
t749 = 0.1e1 / (((t847 * t948 + (0.2e1 * t877 - 0.1e1) * t989) * pkin(1) * t850 + t979) * t853 * t1007 + pkin(1) * t961 * t997 - t902 * t881 + (0.2e1 * (t902 / 0.2e1 - t915 + t947) * t844 + t915) * t912);
t994 = m(3) * t749;
t783 = (t848 * t888 + t851 * t891) * pkin(1) + t894;
t993 = m(3) * t783;
t784 = (t849 * t889 + t852 * t892) * pkin(1) + t894;
t992 = m(3) * t784;
t785 = (t850 * t890 + t853 * t893) * pkin(1) + t894;
t991 = m(3) * t785;
t990 = m(3) * t894;
t857 = cos(t885);
t986 = pkin(2) * t857;
t858 = cos(t886);
t985 = pkin(2) * t858;
t859 = cos(t887);
t984 = pkin(2) * t859;
t983 = pkin(2) * t891;
t982 = pkin(2) * t892;
t981 = pkin(2) * t893;
t837 = t854 * pkin(2);
t838 = t855 * pkin(2);
t839 = t856 * pkin(2);
t980 = qJ(3,1) * t859;
t978 = qJ(3,2) * t858;
t976 = qJ(3,3) * t857;
t803 = t834 + t986;
t974 = t803 * t888;
t807 = t835 + t985;
t973 = t807 * t889;
t811 = t836 + t984;
t972 = t811 * t890;
t971 = t834 * t851;
t970 = t835 * t852;
t969 = t836 * t853;
t968 = t845 * t857;
t967 = t846 * t858;
t966 = t847 * t859;
t964 = t848 * t851;
t962 = t849 * t852;
t960 = t850 * t853;
t959 = t854 * t848;
t958 = t855 * t849;
t957 = t856 * t850;
t956 = t857 * t851;
t955 = t858 * t852;
t954 = t859 * t853;
t883 = t911 / 0.2e1;
t946 = 0.1e1 / 0.2e1 + t883;
t945 = pkin(2) * t980;
t944 = pkin(2) * t978;
t943 = pkin(2) * t976;
t942 = t857 * t965;
t941 = t858 * t963;
t940 = t859 * t961;
t939 = t912 * t950;
t938 = t912 * t949;
t937 = t912 * t948;
t936 = Icges(3,2) + Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t935 = -(2 * t900) - t953;
t934 = -(2 * t901) - t952;
t933 = -(2 * t902) - t951;
t932 = t857 * t939;
t931 = t858 * t938;
t930 = t859 * t937;
t929 = (rSges(3,3) ^ 2) + t911 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t927 = qJ(3,3) * t1006 + t900 + t929;
t771 = t927 * m(3) + t936;
t926 = qJ(3,2) * t1006 + t901 + t929;
t772 = t926 * m(3) + t936;
t925 = qJ(3,1) * t1006 + t902 + t929;
t773 = t925 * m(3) + t936;
t928 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t912 + t936;
t804 = t834 - t986;
t805 = t837 + t976;
t924 = -t804 * t875 - t805 * t950;
t923 = -t804 * t950 + t805 * t875;
t808 = t835 - t985;
t809 = t838 + t978;
t922 = -t808 * t876 - t809 * t949;
t921 = -t808 * t949 + t809 * t876;
t812 = t836 - t984;
t813 = t839 + t980;
t920 = -t812 * t877 - t813 * t948;
t919 = -t812 * t948 + t813 * t877;
t895 = m(2) * rSges(2,2);
t821 = (-qJ(3,3) - rSges(3,3)) * m(3) + t895;
t830 = m(2) * rSges(2,1) + t990;
t918 = (-(t821 * t891 - t830 * t888) * t848 + (t821 * t888 + t830 * t891) * t851) * pkin(1);
t822 = (-qJ(3,2) - rSges(3,3)) * m(3) + t895;
t917 = (-(t822 * t892 - t830 * t889) * t849 + (t822 * t889 + t830 * t892) * t852) * pkin(1);
t823 = (-qJ(3,1) - rSges(3,3)) * m(3) + t895;
t916 = (-(t823 * t893 - t830 * t890) * t850 + (t823 * t890 + t830 * t893) * t853) * pkin(1);
t908 = koppelP(1,1);
t907 = koppelP(2,1);
t906 = koppelP(3,1);
t905 = koppelP(1,2);
t904 = koppelP(2,2);
t903 = koppelP(3,2);
t899 = rSges(4,1);
t898 = rSges(4,2);
t897 = xP(3);
t882 = 0.1e1 + t912;
t874 = -t902 / 0.2e1;
t873 = -t901 / 0.2e1;
t872 = -t900 / 0.2e1;
t867 = cos(t897);
t866 = sin(t897);
t829 = -0.2e1 * t945;
t828 = -0.2e1 * t944;
t827 = -0.2e1 * t943;
t826 = 0.2e1 * t833;
t825 = 0.2e1 * t832;
t824 = 0.2e1 * t831;
t820 = -t979 + t981;
t819 = -t977 + t982;
t818 = -t975 + t983;
t817 = t856 * t937;
t816 = t855 * t938;
t815 = t854 * t939;
t814 = t839 - t980;
t810 = t838 - t978;
t806 = t837 - t976;
t802 = -t866 * t905 + t867 * t908;
t801 = -t866 * t904 + t867 * t907;
t800 = -t866 * t903 + t867 * t906;
t799 = -t866 * t908 - t867 * t905;
t798 = -t866 * t907 - t867 * t904;
t797 = -t866 * t906 - t867 * t903;
t793 = m(4) * (-t866 * t898 + t867 * t899);
t792 = m(4) * (-t866 * t899 - t867 * t898);
t791 = t856 * t862 + t829;
t790 = t855 * t861 + t828;
t789 = t854 * t860 + t827;
t788 = t847 * t856 + t829;
t787 = t846 * t855 + t828;
t786 = t845 * t854 + t827;
t782 = t788 * t893;
t781 = t787 * t892;
t780 = t786 * t891;
t770 = t814 * t893 - t972;
t769 = t811 * t893 + t890 * t814;
t768 = t810 * t892 - t973;
t767 = t807 * t892 + t889 * t810;
t766 = t806 * t891 - t974;
t765 = t803 * t891 + t888 * t806;
t764 = (t859 * t862 + t826) * t893 + t791 * t890;
t763 = (t858 * t861 + t825) * t892 + t790 * t889;
t762 = (t857 * t860 + t824) * t891 + t789 * t888;
t761 = (t826 + t966) * t893 + t890 * t788;
t760 = (t825 + t967) * t892 + t889 * t787;
t759 = (t824 + t968) * t891 + t888 * t786;
t758 = t791 * t893 + ((t883 + t874) * t859 + t833) * t1000;
t757 = t790 * t892 + ((t883 + t873) * t858 + t832) * t1001;
t756 = t789 * t891 + ((t883 + t872) * t857 + t831) * t1002;
t755 = t916 + t773;
t754 = t917 + t772;
t753 = t918 + t771;
t752 = (t912 + t925) * m(3) + 0.2e1 * t916 + t928;
t751 = (t912 + t926) * m(3) + 0.2e1 * t917 + t928;
t750 = (t912 + t927) * m(3) + 0.2e1 * t918 + t928;
t746 = (-t954 + t957) * qJ(3,1) + (-t782 * t844 - t761 * t960 + (t856 * t911 + t856 - t945) * t893 + (-(t1003 - t966) * t844 - qJ(3,1) * t814) * t890) * pkin(1);
t745 = (-t955 + t958) * qJ(3,2) + (-t781 * t843 - t760 * t962 + (t855 * t911 + t855 - t944) * t892 + (-(t1004 - t967) * t843 - qJ(3,2) * t810) * t889) * pkin(1);
t744 = (-t956 + t959) * qJ(3,3) + (-t780 * t842 - t759 * t964 + (t854 * t911 + t854 - t943) * t891 + (-(t1005 - t968) * t842 - qJ(3,3) * t806) * t888) * pkin(1);
t743 = -t940 - t969 + (t761 * t844 - (t782 + ((t874 + t946) * t859 + t833) * t1000) * t960 - (t859 * t911 + t833 + t859) * t893 + qJ(3,1) * t972) * pkin(1);
t742 = -t941 - t970 + (t760 * t843 - (t781 + ((t873 + t946) * t858 + t832) * t1001) * t962 - (t858 * t911 + t832 + t858) * t892 + qJ(3,2) * t973) * pkin(1);
t741 = -t942 - t971 + (t759 * t842 - (t780 + ((t872 + t946) * t857 + t831) * t1002) * t964 - (t857 * t911 + t831 + t857) * t891 + qJ(3,3) * t974) * pkin(1);
t740 = (-t1008 * t856 + t829 + t930) * t853 + (t933 * t859 + t817 + t826) * t850 + (-t770 * t844 - t769 * t960 + (t839 - 0.2e1 * t980) * t893 + t890 * t836) * pkin(1);
t739 = (-t1009 * t855 + t828 + t931) * t852 + (t934 * t858 + t816 + t825) * t849 + (-t768 * t843 - t767 * t962 + (t838 - 0.2e1 * t978) * t892 + t889 * t835) * pkin(1);
t738 = (-t1010 * t854 + t827 + t932) * t851 + (t935 * t857 + t815 + t824) * t848 + (-t766 * t842 - t765 * t964 + (t837 - 0.2e1 * t976) * t891 + t888 * t834) * pkin(1);
t737 = (t1008 * t859 + t1003 + t817) * t853 + (t933 * t856 + t829 - t930) * t850 + (-t770 * t960 + t769 * t844 + t836 * t997 + (-t979 - t981) * t859) * pkin(1);
t736 = (t1009 * t858 + t1004 + t816) * t852 + (t934 * t855 + t828 - t931) * t849 + (-t768 * t962 + t767 * t843 + t835 * t998 + (-t977 - t982) * t858) * pkin(1);
t735 = (t1010 * t857 + t1005 + t815) * t851 + (t935 * t854 + t827 - t932) * t848 + (-t766 * t964 + t765 * t842 + t834 * t999 + (-t975 - t983) * t857) * pkin(1);
t734 = -t882 * t940 - t969 + (t758 * t960 - t764 * t844 + t820 * t811) * pkin(1) + ((t920 - t984) * t853 + t919 * t850) * t912;
t733 = -t882 * t941 - t970 + (t757 * t962 - t763 * t843 + t807 * t819) * pkin(1) + ((t922 - t985) * t852 + t921 * t849) * t912;
t732 = -t882 * t942 - t971 + (t756 * t964 - t762 * t842 + t803 * t818) * pkin(1) + ((t924 - t986) * t851 + t923 * t848) * t912;
t731 = (t882 * t957 - t954) * qJ(3,1) + (t758 * t844 + t764 * t960 - t814 * t820) * pkin(1) + ((-t919 + t839) * t853 + t920 * t850) * t912;
t730 = (t882 * t958 - t955) * qJ(3,2) + (t757 * t843 + t763 * t962 - t810 * t819) * pkin(1) + ((-t921 + t838) * t852 + t922 * t849) * t912;
t729 = (t882 * t959 - t956) * qJ(3,3) + (t756 * t842 + t762 * t964 - t806 * t818) * pkin(1) + ((-t923 + t837) * t851 + t924 * t848) * t912;
t728 = (t743 * t802 + t746 * t799) * t749;
t727 = (t742 * t801 + t745 * t798) * t748;
t726 = (t741 * t800 + t744 * t797) * t747;
t725 = (t737 * t802 + t740 * t799) * t749;
t724 = (t736 * t801 + t739 * t798) * t748;
t723 = (t735 * t800 + t738 * t797) * t747;
t722 = (t731 * t799 + t734 * t802) * t749;
t721 = (t730 * t798 + t733 * t801) * t748;
t720 = (t729 * t797 + t732 * t800) * t747;
t719 = (-t731 * t894 - t746 * t785 + t740) * t994;
t718 = (-t730 * t894 - t745 * t784 + t739) * t995;
t717 = (-t729 * t894 - t744 * t783 + t738) * t996;
t716 = (-t734 * t894 - t743 * t785 + t737) * t994;
t715 = (-t733 * t894 - t742 * t784 + t736) * t995;
t714 = (-t732 * t894 - t741 * t783 + t735) * t996;
t713 = (t731 * t773 - t740 * t990 + t746 * t755) * t749;
t712 = (t730 * t772 - t739 * t990 + t745 * t754) * t748;
t711 = (t729 * t771 - t738 * t990 + t744 * t753) * t747;
t710 = (t734 * t773 - t737 * t990 + t743 * t755) * t749;
t709 = (t733 * t772 - t736 * t990 + t742 * t754) * t748;
t708 = (t732 * t771 - t735 * t990 + t741 * t753) * t747;
t707 = (t731 * t755 - t740 * t991 + t746 * t752) * t749;
t706 = (t730 * t754 - t739 * t992 + t745 * t751) * t748;
t705 = (t729 * t753 - t738 * t993 + t744 * t750) * t747;
t704 = (t734 * t755 - t737 * t991 + t743 * t752) * t749;
t703 = (t733 * t754 - t736 * t992 + t742 * t751) * t748;
t702 = (t732 * t753 - t735 * t993 + t741 * t750) * t747;
t701 = (-t722 * t894 - t728 * t785 + t725) * m(3);
t700 = (-t721 * t894 - t727 * t784 + t724) * m(3);
t699 = (-t720 * t894 - t726 * t783 + t723) * m(3);
t698 = t722 * t773 - t725 * t990 + t728 * t755;
t697 = t721 * t772 - t724 * t990 + t727 * t754;
t696 = t720 * t771 - t723 * t990 + t726 * t753;
t695 = t722 * t755 - t725 * t991 + t728 * t752;
t694 = t721 * t754 - t724 * t992 + t727 * t751;
t693 = t720 * t753 - t723 * t993 + t726 * t750;
t1 = [m(4) + (t707 * t746 + t713 * t731 + t719 * t740) * t749 + (t706 * t745 + t712 * t730 + t718 * t739) * t748 + (t705 * t744 + t711 * t729 + t717 * t738) * t747 (t707 * t743 + t713 * t734 + t719 * t737) * t749 + (t706 * t742 + t712 * t733 + t718 * t736) * t748 + (t705 * t741 + t711 * t732 + t717 * t735) * t747, t705 * t726 + t706 * t727 + t707 * t728 + t711 * t720 + t712 * t721 + t713 * t722 + t717 * t723 + t718 * t724 + t719 * t725 + t792; (t704 * t746 + t710 * t731 + t716 * t740) * t749 + (t703 * t745 + t709 * t730 + t715 * t739) * t748 + (t702 * t744 + t708 * t729 + t714 * t738) * t747, m(4) + (t704 * t743 + t710 * t734 + t716 * t737) * t749 + (t703 * t742 + t709 * t733 + t715 * t736) * t748 + (t702 * t741 + t708 * t732 + t714 * t735) * t747, t702 * t726 + t703 * t727 + t704 * t728 + t708 * t720 + t709 * t721 + t710 * t722 + t714 * t723 + t715 * t724 + t716 * t725 + t793; t792 + (t695 * t746 + t698 * t731 + t701 * t740) * t749 + (t694 * t745 + t697 * t730 + t700 * t739) * t748 + (t693 * t744 + t696 * t729 + t699 * t738) * t747, t793 + (t695 * t743 + t698 * t734 + t701 * t737) * t749 + (t694 * t742 + t697 * t733 + t700 * t736) * t748 + (t693 * t741 + t696 * t732 + t699 * t735) * t747, t695 * t728 + t698 * t722 + t701 * t725 + t694 * t727 + t697 * t721 + t700 * t724 + t693 * t726 + t696 * t720 + t699 * t723 + Icges(4,3) + m(4) * (t898 ^ 2 + t899 ^ 2);];
MX  = t1;
