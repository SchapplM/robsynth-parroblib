% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
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
%   see P3PRR1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PRR1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(8,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_inertia_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_inertia_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'P3PRR1A0_inertia_para_pf_mdp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:41
% EndTime: 2019-05-03 14:47:42
% DurationCPUTime: 0.82s
% Computational Cost: add. (495->122), mult. (1087->245), div. (216->8), fcn. (1190->14), ass. (0->95)
t610 = sin(qJ(2,1));
t603 = 0.1e1 / t610;
t614 = xP(3);
t597 = sin(t614);
t598 = cos(t614);
t617 = koppelP(1,2);
t620 = koppelP(1,1);
t586 = t597 * t620 + t598 * t617;
t607 = legFrame(1,3);
t593 = sin(t607);
t596 = cos(t607);
t623 = t597 * t617 - t598 * t620;
t670 = t586 * t596 + t623 * t593;
t661 = t670 * t603;
t609 = sin(qJ(2,2));
t601 = 0.1e1 / t609;
t616 = koppelP(2,2);
t619 = koppelP(2,1);
t585 = t597 * t619 + t598 * t616;
t606 = legFrame(2,3);
t592 = sin(t606);
t595 = cos(t606);
t624 = t597 * t616 - t598 * t619;
t671 = t585 * t595 + t624 * t592;
t659 = t671 * t601;
t608 = sin(qJ(2,3));
t599 = 0.1e1 / t608;
t615 = koppelP(3,2);
t618 = koppelP(3,1);
t584 = t597 * t618 + t598 * t615;
t605 = legFrame(3,3);
t591 = sin(t605);
t594 = cos(t605);
t625 = t597 * t615 - t598 * t618;
t672 = t584 * t594 + t625 * t591;
t660 = t672 * t599;
t621 = 1 / pkin(2);
t669 = 2 * t621;
t668 = MDP(2) / pkin(2) ^ 2;
t611 = cos(qJ(2,3));
t566 = (t584 * t591 - t594 * t625) * t608 - t672 * t611;
t572 = t621 * t660;
t667 = t566 * t572;
t612 = cos(qJ(2,2));
t567 = (t585 * t592 - t595 * t624) * t609 - t671 * t612;
t573 = t621 * t659;
t666 = t567 * t573;
t613 = cos(qJ(2,1));
t568 = (t586 * t593 - t596 * t623) * t610 - t670 * t613;
t574 = t621 * t661;
t665 = t568 * t574;
t578 = -t591 * t608 + t594 * t611;
t579 = t591 * t611 + t594 * t608;
t569 = (-t578 * t584 - t579 * t625) * t599;
t664 = t569 * t672;
t580 = -t592 * t609 + t595 * t612;
t581 = t592 * t612 + t595 * t609;
t570 = (-t580 * t585 - t581 * t624) * t601;
t663 = t570 * t671;
t582 = -t593 * t610 + t596 * t613;
t583 = t593 * t613 + t596 * t610;
t571 = (-t582 * t586 - t583 * t623) * t603;
t662 = t571 * t670;
t600 = 0.1e1 / t608 ^ 2;
t658 = t579 * t600;
t602 = 0.1e1 / t609 ^ 2;
t657 = t581 * t602;
t604 = 0.1e1 / t610 ^ 2;
t656 = t583 * t604;
t652 = t591 * t599;
t651 = t592 * t601;
t650 = t593 * t603;
t649 = t594 * t599;
t648 = t594 * t600;
t647 = t595 * t601;
t646 = t595 * t602;
t645 = t596 * t603;
t644 = t596 * t604;
t643 = t599 * t611;
t642 = t600 * t611;
t641 = t601 * t612;
t640 = t602 * t612;
t639 = t603 * t613;
t638 = t604 * t613;
t637 = t670 * t638;
t636 = t672 * t642;
t635 = t671 * t640;
t634 = t591 * t642;
t633 = t592 * t640;
t632 = t593 * t638;
t631 = t594 * t642;
t630 = t595 * t640;
t629 = t596 * t638;
t590 = (t597 ^ 2 + t598 ^ 2) * MDP(8);
t1 = [(t578 ^ 2 * t600 + t580 ^ 2 * t602 + t582 ^ 2 * t604) * MDP(1) + t590 + (t594 ^ 2 * t600 + t595 ^ 2 * t602 + t596 ^ 2 * t604) * t668 + ((-t578 * t631 - t580 * t630 - t582 * t629) * MDP(3) + (t578 * t649 + t580 * t647 + t582 * t645) * MDP(4)) * t669; (t578 * t658 + t580 * t657 + t582 * t656) * MDP(1) + (t591 * t648 + t592 * t646 + t593 * t644) * t668 + ((-t578 * t634 - t579 * t631 - t580 * t633 - t581 * t630 - t582 * t632 - t583 * t629) * MDP(3) + (t578 * t652 + t579 * t649 + t580 * t651 + t581 * t647 + t582 * t650 + t583 * t645) * MDP(4)) * t621; (t579 ^ 2 * t600 + t581 ^ 2 * t602 + t583 ^ 2 * t604) * MDP(1) + t590 + (t591 ^ 2 * t600 + t592 ^ 2 * t602 + t593 ^ 2 * t604) * t668 + ((-t579 * t634 - t581 * t633 - t583 * t632) * MDP(3) + (t579 * t652 + t581 * t651 + t583 * t650) * MDP(4)) * t669; (t566 * t578 * t600 + t567 * t580 * t602 + t568 * t582 * t604) * MDP(1) - t597 * MDP(6) - t598 * MDP(7) + (-t644 * t670 - t646 * t671 - t648 * t672) * t668 + ((-t566 * t631 - t567 * t630 - t568 * t629 + t578 * t636 + t580 * t635 + t582 * t637) * MDP(3) + (t566 * t649 + t567 * t647 + t568 * t645 - t578 * t660 - t580 * t659 - t582 * t661) * MDP(4)) * t621; (t566 * t658 + t567 * t657 + t568 * t656) * MDP(1) + t598 * MDP(6) - t597 * MDP(7) + (-t591 * t600 * t672 - t592 * t602 * t671 - t593 * t604 * t670) * t668 + ((-t566 * t634 - t567 * t633 - t568 * t632 + t579 * t636 + t581 * t635 + t583 * t637) * MDP(3) + (t566 * t652 + t567 * t651 + t568 * t650 - t579 * t660 - t581 * t659 - t583 * t661) * MDP(4)) * t621; (t566 * t569 * t599 + t567 * t570 * t601 + t568 * t571 * t603) * MDP(1) + (t639 * t665 + t641 * t666 + t643 * t667) * MDP(3) + (-t665 - t666 - t667) * MDP(4) + MDP(5) + ((t572 * t660 + t573 * t659 + t574 * t661) * MDP(2) + (t639 * t662 + t641 * t663 + t643 * t664) * MDP(3) + (-t662 - t663 - t664) * MDP(4)) * t621;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;
