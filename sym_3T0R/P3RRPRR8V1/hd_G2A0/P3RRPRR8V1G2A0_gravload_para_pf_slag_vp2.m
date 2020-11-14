% Calculate Gravitation load for parallel robot
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:19
% EndTime: 2020-08-06 19:59:19
% DurationCPUTime: 0.37s
% Computational Cost: add. (555->115), mult. (789->190), div. (30->9), fcn. (600->23), ass. (0->86)
t811 = sin(pkin(5));
t868 = pkin(2) * t811;
t816 = legFrame(3,2);
t802 = sin(t816);
t805 = cos(t816);
t788 = t805 * g(1) - t802 * g(2);
t813 = pkin(4) + qJ(3,3);
t808 = 0.1e1 / t813;
t820 = sin(qJ(1,3));
t826 = cos(qJ(1,3));
t867 = (t820 * g(3) - t826 * t788) * t808;
t817 = legFrame(2,2);
t803 = sin(t817);
t806 = cos(t817);
t789 = t806 * g(1) - t803 * g(2);
t814 = pkin(4) + qJ(3,2);
t809 = 0.1e1 / t814;
t822 = sin(qJ(1,2));
t828 = cos(qJ(1,2));
t866 = (t822 * g(3) - t828 * t789) * t809;
t818 = legFrame(1,2);
t804 = sin(t818);
t807 = cos(t818);
t790 = t807 * g(1) - t804 * g(2);
t815 = pkin(4) + qJ(3,1);
t810 = 0.1e1 / t815;
t824 = sin(qJ(1,1));
t830 = cos(qJ(1,1));
t865 = (t824 * g(3) - t830 * t790) * t810;
t812 = cos(pkin(5));
t784 = m(3) * pkin(1) + mrSges(3,1) * t812 - mrSges(3,2) * t811 + mrSges(2,1);
t785 = t802 * g(1) + t805 * g(2);
t794 = t811 * mrSges(3,1) + t812 * mrSges(3,2) + mrSges(2,2);
t819 = sin(qJ(2,3));
t825 = cos(qJ(2,3));
t839 = g(3) * t826 + t788 * t820;
t842 = t825 * pkin(1) + pkin(2) * cos(qJ(2,3) + pkin(5));
t864 = 0.1e1 / t842 * ((t839 * t784 + t785 * t794) * t819 - t825 * (t785 * t784 - t839 * t794));
t786 = t803 * g(1) + t806 * g(2);
t821 = sin(qJ(2,2));
t827 = cos(qJ(2,2));
t838 = g(3) * t828 + t789 * t822;
t841 = t827 * pkin(1) + pkin(2) * cos(qJ(2,2) + pkin(5));
t863 = 0.1e1 / t841 * ((t838 * t784 + t786 * t794) * t821 - t827 * (t786 * t784 - t838 * t794));
t787 = t804 * g(1) + t807 * g(2);
t823 = sin(qJ(2,1));
t829 = cos(qJ(2,1));
t837 = g(3) * t830 + t790 * t824;
t840 = t829 * pkin(1) + pkin(2) * cos(qJ(2,1) + pkin(5));
t862 = 0.1e1 / t840 * ((t837 * t784 + t787 * t794) * t823 - t829 * (t787 * t784 - t837 * t794));
t798 = t812 * pkin(2) + pkin(1);
t861 = t798 * t805;
t860 = t798 * t806;
t859 = t798 * t807;
t858 = t802 * t798;
t857 = t803 * t798;
t856 = t804 * t798;
t852 = mrSges(2,3) + mrSges(3,3) - mrSges(1,2);
t795 = m(3) * qJ(3,3) + t852;
t833 = -t825 * t784 + t794 * t819 - mrSges(1,1);
t855 = t808 * ((-t795 * t820 + t833 * t826) * t788 + (-t795 * t826 - t833 * t820) * g(3));
t796 = m(3) * qJ(3,2) + t852;
t832 = -t827 * t784 + t794 * t821 - mrSges(1,1);
t854 = t809 * ((-t796 * t822 + t832 * t828) * t789 + (-t796 * t828 - t832 * t822) * g(3));
t797 = m(3) * qJ(3,1) + t852;
t831 = -t829 * t784 + t794 * t823 - mrSges(1,1);
t853 = t810 * ((-t797 * t824 + t831 * t830) * t790 + (-t797 * t830 - t831 * t824) * g(3));
t851 = t805 * t868;
t850 = t806 * t868;
t849 = t807 * t868;
t848 = t802 * t868;
t847 = t803 * t868;
t846 = t804 * t868;
t836 = t798 * t825 - t819 * t868;
t845 = 0.1e1 / t836 * t855;
t835 = t798 * t827 - t821 * t868;
t844 = 0.1e1 / t835 * t854;
t834 = t798 * t829 - t823 * t868;
t843 = 0.1e1 / t834 * t853;
t783 = t823 * t798 + t829 * t868;
t782 = t821 * t798 + t827 * t868;
t781 = t819 * t798 + t825 * t868;
t774 = -t830 * t815 + t834 * t824;
t773 = -t828 * t814 + t835 * t822;
t772 = -t826 * t813 + t836 * t820;
t1 = [((t824 * t859 + t846) * t829 + (-t824 * t849 + t856) * t823) * t843 + t804 * t862 + ((t822 * t860 + t847) * t827 + (-t822 * t850 + t857) * t821) * t844 + t803 * t863 + ((t820 * t861 + t848) * t825 + (-t820 * t851 + t858) * t819) * t845 + t802 * t864 - g(1) * m(4) + (-(t774 * t807 + t783 * t804) * t865 - (t773 * t806 + t782 * t803) * t866 - (t772 * t805 + t781 * t802) * t867) * m(3); ((-t824 * t856 + t849) * t829 + t823 * (t824 * t846 + t859)) * t843 + t807 * t862 + ((-t822 * t857 + t850) * t827 + t821 * (t822 * t847 + t860)) * t844 + t806 * t863 + ((-t820 * t858 + t851) * t825 + t819 * (t820 * t848 + t861)) * t845 + t805 * t864 - g(2) * m(4) + (-(-t774 * t804 + t783 * t807) * t865 - (-t773 * t803 + t782 * t806) * t866 - (-t772 * t802 + t781 * t805) * t867) * m(3); t826 * t855 + t828 * t854 + t830 * t853 - g(3) * m(4) + (-(t824 * t815 + t840 * t830) * t865 - (t822 * t814 + t841 * t828) * t866 - (t820 * t813 + t842 * t826) * t867) * m(3);];
taugX  = t1;
