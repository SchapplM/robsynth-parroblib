% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [8x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PRR1G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRR1G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:42
% EndTime: 2019-05-03 14:47:43
% DurationCPUTime: 0.75s
% Computational Cost: add. (865->109), mult. (2112->213), div. (369->8), fcn. (1612->14), ass. (0->89)
t733 = xP(3);
t709 = sin(t733);
t710 = cos(t733);
t734 = koppelP(3,2);
t737 = koppelP(3,1);
t697 = t709 * t737 + t710 * t734;
t721 = legFrame(3,3);
t703 = sin(t721);
t706 = cos(t721);
t753 = t709 * t734 - t710 * t737;
t689 = t697 * t706 + t703 * t753;
t735 = koppelP(2,2);
t738 = koppelP(2,1);
t698 = t709 * t738 + t710 * t735;
t722 = legFrame(2,3);
t704 = sin(t722);
t707 = cos(t722);
t752 = t709 * t735 - t710 * t738;
t690 = t698 * t707 + t704 * t752;
t736 = koppelP(1,2);
t739 = koppelP(1,1);
t699 = t709 * t739 + t710 * t736;
t723 = legFrame(1,3);
t705 = sin(t723);
t708 = cos(t723);
t751 = t709 * t736 - t710 * t739;
t688 = t699 * t708 + t705 * t751;
t730 = xDP(3);
t731 = xDP(2);
t732 = xDP(1);
t685 = t689 * t730 - t703 * t731 - t732 * t706;
t724 = sin(qJ(2,3));
t775 = t685 ^ 2 / t724 ^ 2;
t686 = t690 * t730 - t704 * t731 - t732 * t707;
t725 = sin(qJ(2,2));
t773 = t686 ^ 2 / t725 ^ 2;
t687 = t688 * t730 - t705 * t731 - t732 * t708;
t726 = sin(qJ(2,1));
t771 = t687 ^ 2 / t726 ^ 2;
t711 = 0.1e1 / t724;
t714 = 0.1e1 / t725;
t717 = 0.1e1 / t726;
t727 = cos(qJ(2,3));
t691 = -t703 * t724 + t706 * t727;
t692 = t703 * t727 + t724 * t706;
t720 = t730 ^ 2;
t740 = 0.1e1 / pkin(2);
t774 = t711 * t775;
t673 = t740 * t774 + (t691 * t753 - t692 * t697) * t720 * t711;
t781 = t673 * t711;
t728 = cos(qJ(2,2));
t693 = -t704 * t725 + t707 * t728;
t694 = t704 * t728 + t725 * t707;
t772 = t714 * t773;
t674 = t740 * t772 + (t693 * t752 - t694 * t698) * t720 * t714;
t780 = t674 * t714;
t729 = cos(qJ(2,1));
t695 = -t705 * t726 + t708 * t729;
t696 = t705 * t729 + t726 * t708;
t770 = t717 * t771;
t675 = t740 * t770 + (t695 * t751 - t696 * t699) * t720 * t717;
t779 = t675 * t717;
t741 = 0.1e1 / pkin(2) ^ 2;
t756 = t703 * t697 - t706 * t753;
t759 = t727 * t774;
t766 = t720 * t740;
t676 = t711 * t756 * t766 - t741 * t759;
t778 = t676 * t711;
t755 = t704 * t698 - t707 * t752;
t758 = t728 * t772;
t677 = t714 * t755 * t766 - t741 * t758;
t777 = t677 * t714;
t754 = t705 * t699 - t708 * t751;
t757 = t729 * t770;
t678 = t717 * t754 * t766 - t741 * t757;
t776 = t678 * t717;
t769 = t711 * t727;
t768 = t714 * t728;
t767 = t717 * t729;
t765 = t673 * t769;
t764 = t674 * t768;
t763 = t675 * t767;
t762 = t676 * t769;
t761 = t677 * t768;
t760 = t678 * t767;
t681 = -t688 * t729 + t754 * t726;
t680 = -t690 * t728 + t755 * t725;
t679 = -t689 * t727 + t756 * t724;
t1 = [(t691 * t781 + t693 * t780 + t695 * t779) * MDP(1) + (t691 * t762 + t693 * t761 + t695 * t760) * MDP(3) + (-t691 * t676 - t693 * t677 - t695 * t678) * MDP(4) + (-t710 * MDP(6) + t709 * MDP(7)) * t720 + ((-t691 * t775 - t693 * t773 - t695 * t771) * MDP(3) + (-t691 * t759 - t693 * t758 - t695 * t757) * MDP(4)) * t741 + ((-t706 * t778 - t707 * t777 - t708 * t776) * MDP(2) + (-t706 * t765 - t707 * t764 - t708 * t763) * MDP(3) + (t673 * t706 + t674 * t707 + t675 * t708) * MDP(4)) * t740; (t692 * t781 + t694 * t780 + t696 * t779) * MDP(1) + (t692 * t762 + t694 * t761 + t696 * t760) * MDP(3) + (-t692 * t676 - t694 * t677 - t696 * t678) * MDP(4) + (-t709 * MDP(6) - t710 * MDP(7)) * t720 + ((-t692 * t775 - t694 * t773 - t696 * t771) * MDP(3) + (-t692 * t759 - t694 * t758 - t696 * t757) * MDP(4)) * t741 + ((-t703 * t778 - t704 * t777 - t705 * t776) * MDP(2) + (-t703 * t765 - t704 * t764 - t705 * t763) * MDP(3) + (t673 * t703 + t674 * t704 + t675 * t705) * MDP(4)) * t740; (t679 * t781 + t680 * t780 + t681 * t779) * MDP(1) + (t679 * t762 + t680 * t761 + t681 * t760) * MDP(3) + (-t679 * t676 - t680 * t677 - t681 * t678) * MDP(4) + ((-t679 * t775 - t680 * t773 - t681 * t771) * MDP(3) + (-t679 * t759 - t680 * t758 - t681 * t757) * MDP(4)) * t741 + ((t688 * t776 + t689 * t778 + t690 * t777) * MDP(2) + (t688 * t763 + t689 * t765 + t690 * t764) * MDP(3) + (-t673 * t689 - t674 * t690 - t675 * t688) * MDP(4)) * t740;];
taucX  = t1;
