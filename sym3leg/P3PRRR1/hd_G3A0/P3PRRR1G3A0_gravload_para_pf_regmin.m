% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PRRR1G3P3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR1G3P3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:58
% EndTime: 2020-03-09 21:06:59
% DurationCPUTime: 0.28s
% Computational Cost: add. (865->82), mult. (468->156), div. (126->5), fcn. (558->18), ass. (0->78)
t780 = legFrame(3,2);
t771 = sin(t780);
t777 = pkin(7) + qJ(2,3);
t768 = qJ(3,3) + t777;
t756 = sin(t768);
t759 = cos(t768);
t762 = sin(t777);
t765 = cos(t777);
t805 = 0.1e1 / (t756 * t765 - t762 * t759);
t811 = t771 * t805;
t781 = legFrame(2,2);
t772 = sin(t781);
t778 = pkin(7) + qJ(2,2);
t769 = qJ(3,2) + t778;
t757 = sin(t769);
t760 = cos(t769);
t763 = sin(t778);
t766 = cos(t778);
t804 = 0.1e1 / (t757 * t766 - t763 * t760);
t810 = t772 * t804;
t782 = legFrame(1,2);
t773 = sin(t782);
t779 = pkin(7) + qJ(2,1);
t770 = qJ(3,1) + t779;
t758 = sin(t770);
t761 = cos(t770);
t764 = sin(t779);
t767 = cos(t779);
t803 = 0.1e1 / (t758 * t767 - t764 * t761);
t809 = t773 * t803;
t774 = cos(t780);
t808 = t774 * t805;
t775 = cos(t781);
t807 = t775 * t804;
t776 = cos(t782);
t806 = t776 * t803;
t802 = t805 * (pkin(2) * t762 + pkin(3) * t756);
t801 = t804 * (pkin(2) * t763 + pkin(3) * t757);
t800 = t803 * (pkin(2) * t764 + pkin(3) * t758);
t799 = t805 * t756;
t798 = t804 * t757;
t797 = t803 * t758;
t747 = pkin(2) * t765 + pkin(3) * t759;
t796 = t747 * t811;
t795 = t759 * t808;
t748 = pkin(2) * t766 + pkin(3) * t760;
t794 = t748 * t810;
t793 = t760 * t807;
t749 = pkin(2) * t767 + pkin(3) * t761;
t792 = t749 * t809;
t791 = t761 * t806;
t790 = t747 * t808;
t789 = t759 * t811;
t788 = t748 * t807;
t787 = t760 * t810;
t786 = t749 * t806;
t785 = t761 * t809;
t784 = 0.1e1 / pkin(2);
t783 = 0.1e1 / pkin(3);
t755 = t776 * g(1) - t773 * g(2);
t754 = t775 * g(1) - t772 * g(2);
t753 = t774 * g(1) - t771 * g(2);
t752 = -t773 * g(1) - t776 * g(2);
t751 = -t772 * g(1) - t775 * g(2);
t750 = -t771 * g(1) - t774 * g(2);
t734 = -g(3) * t764 + t755 * t767;
t733 = g(3) * t767 + t755 * t764;
t732 = -g(3) * t763 + t754 * t766;
t731 = g(3) * t766 + t754 * t763;
t730 = -g(3) * t762 + t753 * t765;
t729 = g(3) * t765 + t753 * t762;
t728 = -g(3) * t758 + t755 * t761;
t727 = g(3) * t761 + t755 * t758;
t726 = -g(3) * t757 + t754 * t760;
t725 = g(3) * t760 + t754 * t757;
t724 = -g(3) * t756 + t753 * t759;
t723 = g(3) * t759 + t753 * t756;
t1 = [t771 * t750 + t772 * t751 + t773 * t752, 0, (t729 * t795 + t731 * t793 + t733 * t791) * t784, (t730 * t795 + t732 * t793 + t734 * t791) * t784, 0, (t723 * t795 + t725 * t793 + t727 * t791 + (-t723 * t790 - t725 * t788 - t727 * t786) * t783) * t784, (t724 * t795 + t726 * t793 + t728 * t791 + (-t724 * t790 - t726 * t788 - t728 * t786) * t783) * t784, -g(1); t774 * t750 + t775 * t751 + t776 * t752, 0, (-t729 * t789 - t731 * t787 - t733 * t785) * t784, (-t730 * t789 - t732 * t787 - t734 * t785) * t784, 0, (-t723 * t789 - t725 * t787 - t727 * t785 + (t723 * t796 + t725 * t794 + t727 * t792) * t783) * t784, (-t724 * t789 - t726 * t787 - t728 * t785 + (t724 * t796 + t726 * t794 + t728 * t792) * t783) * t784, -g(2); 0, 0, (-t729 * t799 - t731 * t798 - t733 * t797) * t784, (-t730 * t799 - t732 * t798 - t734 * t797) * t784, 0, (-t723 * t799 - t725 * t798 - t727 * t797 + (t723 * t802 + t725 * t801 + t727 * t800) * t783) * t784, (-t724 * t799 - t726 * t798 - t728 * t797 + (t724 * t802 + t726 * t801 + t728 * t800) * t783) * t784, -g(3);];
tau_reg  = t1;
