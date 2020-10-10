% Calculate vector of centrifugal and coriolis load on the joints for
% P3RRPRR8V2G2A0
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
% Datum: 2020-08-06 21:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:12:04
% EndTime: 2020-08-06 21:12:10
% DurationCPUTime: 6.06s
% Computational Cost: add. (26973->439), mult. (46998->672), div. (4179->12), fcn. (33657->50), ass. (0->316)
t821 = -qJ(3,1) - pkin(5);
t1002 = mrSges(3,2) * t821 + Ifges(3,6);
t757 = mrSges(3,1) * t821 + Ifges(3,5);
t816 = sin(pkin(7));
t817 = cos(pkin(7));
t1017 = t1002 * t816 - t757 * t817;
t820 = -qJ(3,2) - pkin(5);
t1003 = mrSges(3,2) * t820 + Ifges(3,6);
t756 = mrSges(3,1) * t820 + Ifges(3,5);
t1016 = t1003 * t816 - t756 * t817;
t819 = -qJ(3,3) - pkin(5);
t1004 = mrSges(3,2) * t819 + Ifges(3,6);
t755 = mrSges(3,1) * t819 + Ifges(3,5);
t1015 = t1004 * t816 - t755 * t817;
t825 = Ifges(3,2) - Ifges(3,1);
t746 = -pkin(2) * mrSges(3,2) - t816 * t825;
t961 = mrSges(3,1) * t816;
t751 = -pkin(2) * t961 + Ifges(2,4) - Ifges(3,4);
t804 = t817 ^ 2;
t959 = Ifges(3,4) * t804;
t724 = t746 * t817 + t751 + 0.2e1 * t959;
t1014 = -0.2e1 * pkin(2);
t919 = 2 * pkin(3);
t1013 = pkin(1) / 0.2e1;
t1012 = pkin(2) / 0.2e1;
t810 = qJ(2,1) + pkin(7);
t779 = cos(t810);
t764 = pkin(3) * t779;
t841 = cos(qJ(2,1));
t796 = t841 * pkin(2);
t922 = t764 + t796;
t1011 = pkin(1) + t922;
t809 = qJ(2,2) + pkin(7);
t778 = cos(t809);
t763 = pkin(3) * t778;
t839 = cos(qJ(2,2));
t795 = t839 * pkin(2);
t923 = t763 + t795;
t1010 = pkin(1) + t923;
t808 = qJ(2,3) + pkin(7);
t777 = cos(t808);
t762 = pkin(3) * t777;
t837 = cos(qJ(2,3));
t794 = t837 * pkin(2);
t924 = t762 + t794;
t1009 = pkin(1) + t924;
t1001 = 0.2e1 * pkin(1);
t1000 = 0.2e1 * mrSges(3,1);
t999 = -2 * mrSges(3,3);
t738 = 0.1e1 / t924;
t828 = legFrame(3,2);
t784 = sin(t828);
t787 = cos(t828);
t846 = xDP(2);
t847 = xDP(1);
t708 = (t784 * t847 + t787 * t846) * t738;
t705 = t708 ^ 2;
t739 = 0.1e1 / t923;
t829 = legFrame(2,2);
t785 = sin(t829);
t788 = cos(t829);
t709 = (t785 * t847 + t788 * t846) * t739;
t706 = t709 ^ 2;
t740 = 0.1e1 / t922;
t830 = legFrame(1,2);
t786 = sin(t830);
t789 = cos(t830);
t710 = (t786 * t847 + t789 * t846) * t740;
t707 = t710 ^ 2;
t966 = pkin(3) * t816;
t773 = pkin(1) * t966;
t831 = sin(qJ(2,3));
t865 = pkin(3) ^ 2;
t866 = pkin(2) ^ 2;
t965 = pkin(3) * t817;
t772 = pkin(2) * t965;
t990 = 0.2e1 * t772;
t873 = 0.2e1 * t804 * t865 - t865 + t866 + t990;
t721 = t873 * t831 + t773;
t801 = pkin(6) - t819;
t832 = sin(qJ(1,3));
t838 = cos(qJ(1,3));
t891 = pkin(1) * t832 - t801 * t838;
t906 = t832 * t966;
t986 = -0.2e1 * t831;
t728 = t906 * t986 + t891;
t818 = t866 / 0.2e1;
t921 = t772 + t818;
t735 = (t804 - 0.1e1 / 0.2e1) * t865 + t921;
t972 = pkin(1) * t831;
t747 = -t966 + t972;
t765 = pkin(2) + t965;
t909 = t787 * t966;
t942 = t784 * t832;
t945 = t765 * t787;
t948 = t765 * t784;
t903 = pkin(3) * (t804 - 0.1e1);
t936 = t816 * t831;
t969 = pkin(3) * (t832 * t903 + t891 * t936);
t812 = t837 ^ 2;
t989 = 0.2e1 * t812;
t684 = (-t735 * t942 + t765 * t909) * t989 + (t721 * t787 - t728 * t948) * t837 + t784 * t969 + t747 * t945;
t883 = pkin(3) * t936 - t765 * t837;
t731 = 0.1e1 / t883;
t781 = 0.1e1 / t801;
t952 = t731 * t781;
t675 = t684 * t846 * t952;
t912 = t784 * t966;
t939 = t787 * t832;
t685 = (t735 * t939 + t765 * t912) * t989 + (t721 * t784 + t728 * t945) * t837 - t787 * t969 + t747 * t948;
t676 = t685 * t847 * t952;
t725 = t1009 * t838 + t832 * t801;
t845 = xDP(3);
t715 = t725 * t781 * t845;
t669 = -t676 - t675 + t715;
t998 = 0.2e1 * t669;
t833 = sin(qJ(2,2));
t722 = t873 * t833 + t773;
t802 = pkin(6) - t820;
t834 = sin(qJ(1,2));
t840 = cos(qJ(1,2));
t890 = pkin(1) * t834 - t802 * t840;
t905 = t834 * t966;
t985 = -0.2e1 * t833;
t729 = t905 * t985 + t890;
t971 = pkin(1) * t833;
t748 = -t966 + t971;
t908 = t788 * t966;
t941 = t785 * t834;
t944 = t765 * t788;
t947 = t765 * t785;
t935 = t816 * t833;
t968 = pkin(3) * (t834 * t903 + t890 * t935);
t813 = t839 ^ 2;
t988 = 0.2e1 * t813;
t686 = (-t735 * t941 + t765 * t908) * t988 + (t722 * t788 - t729 * t947) * t839 + t785 * t968 + t748 * t944;
t882 = pkin(3) * t935 - t765 * t839;
t732 = 0.1e1 / t882;
t782 = 0.1e1 / t802;
t951 = t732 * t782;
t677 = t686 * t846 * t951;
t911 = t785 * t966;
t938 = t788 * t834;
t687 = (t735 * t938 + t765 * t911) * t988 + (t722 * t785 + t729 * t944) * t839 - t788 * t968 + t748 * t947;
t678 = t687 * t847 * t951;
t726 = t1010 * t840 + t834 * t802;
t716 = t726 * t782 * t845;
t670 = -t678 - t677 + t716;
t997 = 0.2e1 * t670;
t835 = sin(qJ(2,1));
t723 = t873 * t835 + t773;
t803 = pkin(6) - t821;
t836 = sin(qJ(1,1));
t842 = cos(qJ(1,1));
t889 = pkin(1) * t836 - t803 * t842;
t904 = t836 * t966;
t984 = -0.2e1 * t835;
t730 = t904 * t984 + t889;
t970 = pkin(1) * t835;
t749 = -t966 + t970;
t907 = t789 * t966;
t940 = t786 * t836;
t943 = t765 * t789;
t946 = t765 * t786;
t934 = t816 * t835;
t967 = pkin(3) * (t836 * t903 + t889 * t934);
t814 = t841 ^ 2;
t987 = 0.2e1 * t814;
t688 = (-t735 * t940 + t765 * t907) * t987 + (t723 * t789 - t730 * t946) * t841 + t786 * t967 + t749 * t943;
t881 = pkin(3) * t934 - t765 * t841;
t733 = 0.1e1 / t881;
t783 = 0.1e1 / t803;
t950 = t733 * t783;
t679 = t688 * t846 * t950;
t910 = t786 * t966;
t937 = t789 * t836;
t689 = (t735 * t937 + t765 * t910) * t987 + (t723 * t786 + t730 * t943) * t841 - t789 * t967 + t749 * t946;
t680 = t689 * t847 * t950;
t727 = t1011 * t842 + t836 * t803;
t717 = t727 * t783 * t845;
t671 = -t680 - t679 + t717;
t996 = 0.2e1 * t671;
t699 = (-t765 * t942 + t909) * t837 + (t784 * t906 + t945) * t831;
t702 = (t765 * t939 + t912) * t837 + t831 * (-t787 * t906 + t948);
t693 = (t838 * t845 - (t699 * t846 + t702 * t847) * t731) * t781;
t995 = 0.2e1 * t693;
t700 = (-t765 * t941 + t908) * t839 + (t785 * t905 + t944) * t833;
t703 = (t765 * t938 + t911) * t839 + t833 * (-t788 * t905 + t947);
t694 = (t840 * t845 - (t700 * t846 + t703 * t847) * t732) * t782;
t994 = 0.2e1 * t694;
t701 = (-t765 * t940 + t907) * t841 + (t786 * t904 + t943) * t835;
t704 = (t765 * t937 + t910) * t841 + t835 * (-t789 * t904 + t946);
t695 = (t842 * t845 - (t701 * t846 + t704 * t847) * t733) * t783;
t993 = 0.2e1 * t695;
t933 = t825 * t804;
t958 = Ifges(3,4) * t816;
t960 = mrSges(3,2) * t816;
t974 = m(3) * t866;
t714 = 0.2e1 * t933 + (pkin(2) * t1000 + 0.4e1 * t958) * t817 + t974 + t960 * t1014 + Ifges(2,2) - Ifges(2,1) - t825;
t992 = -0.2e1 * t714;
t991 = -0.4e1 * pkin(1) * (t865 / 0.2e1 + t921);
t983 = -0.4e1 * t837;
t982 = -0.4e1 * t839;
t981 = -0.4e1 * t841;
t980 = -4 * pkin(5) - 4 * pkin(6);
t853 = m(3) * pkin(2);
t979 = mrSges(2,2) * pkin(1);
t978 = pkin(5) * mrSges(2,2);
t977 = t738 / 0.2e1;
t976 = t739 / 0.2e1;
t975 = t740 / 0.2e1;
t790 = mrSges(2,1) + t853;
t973 = pkin(1) * t790;
t964 = pkin(3) * t866;
t690 = t695 * pkin(1);
t963 = (mrSges(2,3) + mrSges(3,3)) * pkin(5);
t962 = t865 * pkin(2);
t957 = -t790 * pkin(5) + Ifges(2,5);
t956 = pkin(5) * mrSges(2,1) - Ifges(2,5);
t955 = t724 * t812;
t954 = t724 * t813;
t953 = t724 * t814;
t932 = t831 * t837;
t931 = t833 * t839;
t930 = t835 * t841;
t929 = -0.2e1 * t853;
t928 = pkin(2) * t919;
t691 = pkin(1) * t693;
t664 = t691 + t676 / 0.2e1 + t675 / 0.2e1 - t715 / 0.2e1;
t742 = pkin(2) * t831 + pkin(3) * sin(t808);
t766 = 0.2e1 * t808;
t771 = pkin(2) * t818 + t962;
t860 = 0.2e1 * qJ(2,3);
t807 = t860 + pkin(7);
t858 = 0.2e1 * pkin(7);
t863 = pkin(5) ^ 2;
t867 = pkin(1) ^ 2;
t920 = t865 + t866;
t880 = -(2 * t863) - 0.2e1 * t867 - t920 + ((-4 * pkin(5) - 2 * pkin(6)) * pkin(6));
t895 = -(2 * pkin(3) * t865) - 0.4e1 * t964;
t916 = -0.2e1 * t962;
t918 = -0.2e1 * t964;
t648 = ((cos(qJ(2,3) - pkin(7)) * t918 + cos(t858 + qJ(2,3)) * t916 + t895 * t777 + t771 * t983 + t991) * t705 * t977 + (t1009 * t669 - 0.2e1 * t664 * t762 + t998 * t1013 + (-cos(t766) * t865 - cos(t860) * t866 + t880 + ((t980 - 2 * qJ(3,3)) * qJ(3,3))) * t693 / 0.2e1 + (t664 * t983 + (-cos(t807) - t817) * t693 * t919) * t1012 + (t742 + (sin(t807) * t928 + sin(t766) * t865 + sin(t860) * t866) * t977) * t708 * t801) * t693) * t781;
t750 = t990 + t920;
t660 = (-t705 * t750 / (t794 + (t817 * t837 - t936) * pkin(3)) + (t883 * t693 - t691 + t998) * t693) * t781;
t681 = t693 * t973;
t896 = -m(3) * qJ(3,3) - mrSges(3,3);
t899 = Ifges(2,6) - t978;
t696 = -(t896 * pkin(2) - t1015 + t957) * t831 - t837 * (t1004 * t817 + t755 * t816 + t899);
t888 = mrSges(3,1) * t817 - t960;
t736 = -t853 - t888;
t854 = m(3) * pkin(1);
t887 = mrSges(3,2) * t817 + t961;
t711 = -t736 * t837 - t887 * t831 + t854;
t741 = 0.2e1 * t746;
t745 = (t790 - t960) * t1001;
t852 = m(3) * pkin(5);
t758 = t852 - t896;
t769 = -t978 / 0.2e1 + Ifges(2,6) / 0.2e1;
t859 = -pkin(5) / 0.2e1;
t791 = -qJ(3,3) / 0.2e1 + t859;
t811 = pkin(1) * t1000;
t850 = Ifges(3,6) / 0.2e1;
t851 = Ifges(3,5) / 0.2e1;
t870 = -(m(2) * t863) - (m(2) + m(3)) * t867 - Ifges(2,1) - Ifges(3,2) - Ifges(1,3) + t933;
t871 = 0.2e1 * t724;
t872 = 0.2e1 * mrSges(2,2) + 0.2e1 * t887;
t879 = t888 * t831;
t902 = t705 * t738 * t742;
t915 = 0.4e1 * t959;
t917 = (mrSges(2,2) + t961) * t1001;
t927 = (-t714 * t812 - (t831 * t915 + (t741 * t831 + t811) * t817 + t745) * t837 + t831 * t917 - t819 ^ 2 * m(3) + (qJ(3,3) * t999) + t870) * t660 - t696 * t902 + t711 * t648 + 0.2e1 * (-t751 * t932 - (-mrSges(3,2) * t972 - t958) * t817 - t963) * t660 + t669 * t758 * t995 + (t681 * t986 + (-(t758 * pkin(2) + t1015 + t956) * t837 + ((mrSges(3,2) * t791 + t850) * t817 + (mrSges(3,1) * t791 + t851) * t816 + t769) * t986) * t708 + (t932 * t992 + 0.4e1 * t955 + (-t837 * t872 - 0.2e1 * t879) * pkin(1) - t871) * t693) * t708;
t692 = pkin(1) * t694;
t665 = t692 + t678 / 0.2e1 + t677 / 0.2e1 - t716 / 0.2e1;
t743 = pkin(2) * t833 + pkin(3) * sin(t809);
t767 = 0.2e1 * t809;
t861 = 0.2e1 * qJ(2,2);
t805 = pkin(7) + t861;
t649 = ((cos(qJ(2,2) - pkin(7)) * t918 + cos(t858 + qJ(2,2)) * t916 + t895 * t778 + t771 * t982 + t991) * t706 * t976 + (t1010 * t670 - 0.2e1 * t665 * t763 + t997 * t1013 + (-cos(t767) * t865 - cos(t861) * t866 + t880 + ((t980 - 2 * qJ(3,2)) * qJ(3,2))) * t694 / 0.2e1 + (t665 * t982 + (-cos(t805) - t817) * t694 * t919) * t1012 + (t743 + (sin(t805) * t928 + sin(t767) * t865 + sin(t861) * t866) * t976) * t709 * t802) * t694) * t782;
t661 = (-t706 * t750 / (t795 + (t817 * t839 - t935) * pkin(3)) + (t882 * t694 - t692 + t997) * t694) * t782;
t682 = t694 * t973;
t897 = -m(3) * qJ(3,2) - mrSges(3,3);
t697 = -(t897 * pkin(2) - t1016 + t957) * t833 - t839 * (t1003 * t817 + t756 * t816 + t899);
t712 = -t736 * t839 - t887 * t833 + t854;
t759 = t852 - t897;
t792 = -qJ(3,2) / 0.2e1 + t859;
t878 = t888 * t833;
t901 = t706 * t739 * t743;
t926 = (-t714 * t813 - (t833 * t915 + (t741 * t833 + t811) * t817 + t745) * t839 + t833 * t917 - t820 ^ 2 * m(3) + (qJ(3,2) * t999) + t870) * t661 - t697 * t901 + t712 * t649 + 0.2e1 * (-t751 * t931 - (-mrSges(3,2) * t971 - t958) * t817 - t963) * t661 + t670 * t759 * t994 + (t682 * t985 + (-(t759 * pkin(2) + t1016 + t956) * t839 + ((mrSges(3,2) * t792 + t850) * t817 + (mrSges(3,1) * t792 + t851) * t816 + t769) * t985) * t709 + (t931 * t992 + 0.4e1 * t954 + (-t839 * t872 - 0.2e1 * t878) * pkin(1) - t871) * t694) * t709;
t663 = t690 + t680 / 0.2e1 + t679 / 0.2e1 - t717 / 0.2e1;
t744 = pkin(2) * t835 + pkin(3) * sin(t810);
t768 = 0.2e1 * t810;
t862 = 0.2e1 * qJ(2,1);
t806 = pkin(7) + t862;
t650 = ((cos(qJ(2,1) - pkin(7)) * t918 + cos(qJ(2,1) + t858) * t916 + t895 * t779 + t771 * t981 + t991) * t707 * t975 + (t1011 * t671 - 0.2e1 * t663 * t764 + t996 * t1013 + (-t865 * cos(t768) - t866 * cos(t862) + t880 + ((t980 - 2 * qJ(3,1)) * qJ(3,1))) * t695 / 0.2e1 + (t663 * t981 + (-cos(t806) - t817) * t695 * t919) * t1012 + (t744 + (sin(t806) * t928 + sin(t768) * t865 + sin(t862) * t866) * t975) * t710 * t803) * t695) * t783;
t662 = (-t707 * t750 / (t796 + (t817 * t841 - t934) * pkin(3)) + (t881 * t695 - t690 + t996) * t695) * t783;
t683 = t790 * t690;
t898 = -m(3) * qJ(3,1) - mrSges(3,3);
t698 = -(t898 * pkin(2) - t1017 + t957) * t835 - t841 * (t1002 * t817 + t757 * t816 + t899);
t713 = -t736 * t841 - t887 * t835 + t854;
t760 = t852 - t898;
t793 = -qJ(3,1) / 0.2e1 + t859;
t877 = t888 * t835;
t900 = t707 * t740 * t744;
t925 = (-t714 * t814 - (t835 * t915 + (t741 * t835 + t811) * t817 + t745) * t841 + t835 * t917 - t821 ^ 2 * m(3) + (qJ(3,1) * t999) + t870) * t662 - t698 * t900 + t713 * t650 + 0.2e1 * (-t751 * t930 - (-mrSges(3,2) * t970 - t958) * t817 - t963) * t662 + t671 * t760 * t993 + (t683 * t984 + (-(t760 * pkin(2) + t1017 + t956) * t841 + ((mrSges(3,2) * t793 + t850) * t817 + (mrSges(3,1) * t793 + t851) * t816 + t769) * t984) * t710 + (t930 * t992 + 0.4e1 * t953 + (-t841 * t872 - 0.2e1 * t877) * pkin(1) - t871) * t695) * t710;
t876 = t887 * t837;
t894 = -m(3) * t648 + t660 * t711 + (-t693 * (-m(3) * t819 + mrSges(3,3)) / 0.2e1 + (-t736 * t831 + t876) * t708) * t995;
t875 = t887 * t839;
t893 = -m(3) * t649 + t661 * t712 + (-t694 * (-m(3) * t820 + mrSges(3,3)) / 0.2e1 + (-t736 * t833 + t875) * t709) * t994;
t874 = t887 * t841;
t892 = -m(3) * t650 + t662 * t713 + (-t695 * (-m(3) * t821 + mrSges(3,3)) / 0.2e1 + (-t736 * t835 + t874) * t710) * t993;
t734 = t888 * t1014 - Ifges(2,3) - Ifges(3,3) - t974;
t886 = t738 * (((t669 * t929 + t681) * t831 + (t876 + t879) * (t691 + 0.2e1 * t676 + 0.2e1 * t675 - 0.2e1 * t715) + (-0.2e1 * t955 + (t714 * t831 + t979) * t837 + t724) * t693) * t693 + t660 * t696 - t734 * t902);
t885 = t739 * (((t670 * t929 + t682) * t833 + (t875 + t878) * (t692 + 0.2e1 * t678 + 0.2e1 * t677 - 0.2e1 * t716) + (-0.2e1 * t954 + (t714 * t833 + t979) * t839 + t724) * t694) * t694 + t661 * t697 - t734 * t901);
t884 = t740 * (((t671 * t929 + t683) * t835 + (t874 + t877) * (t690 + 0.2e1 * t680 + 0.2e1 * t679 - 0.2e1 * t717) + (-0.2e1 * t953 + (t714 * t835 + t979) * t841 + t724) * t695) * t695 + t662 * t698 - t734 * t900);
t1 = [t786 * t884 + t785 * t885 + t784 * t886 - (t892 * t689 + t925 * t704) * t950 - (t893 * t687 + t926 * t703) * t951 - (t894 * t685 + t927 * t702) * t952; t789 * t884 + t788 * t885 + t787 * t886 - (t892 * t688 + t925 * t701) * t950 - (t893 * t686 + t926 * t700) * t951 - (t894 * t684 + t927 * t699) * t952; (t892 * t727 + t925 * t842) * t783 + (t893 * t726 + t926 * t840) * t782 + (t894 * t725 + t927 * t838) * t781;];
taucX  = t1;
