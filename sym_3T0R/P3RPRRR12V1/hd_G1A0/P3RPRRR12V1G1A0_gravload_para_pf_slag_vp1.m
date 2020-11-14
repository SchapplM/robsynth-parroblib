% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:21
% EndTime: 2020-08-06 18:21:21
% DurationCPUTime: 0.37s
% Computational Cost: add. (339->91), mult. (492->151), div. (21->7), fcn. (321->18), ass. (0->62)
t642 = pkin(1) + pkin(5);
t641 = rSges(3,2) * g(3);
t587 = (rSges(3,3) + t642) * m(3) + (pkin(1) - rSges(2,2)) * m(2) + m(1) * rSges(1,1);
t613 = legFrame(3,3);
t603 = sin(t613);
t606 = cos(t613);
t588 = -t603 * g(1) + t606 * g(2);
t591 = t606 * g(1) + t603 * g(2);
t629 = m(1) * rSges(1,2);
t594 = -qJ(2,3) * m(3) + (-rSges(2,3) - qJ(2,3)) * m(2) + t629;
t616 = sin(qJ(3,3));
t600 = t616 * pkin(3) + qJ(2,3);
t597 = 0.1e1 / t600;
t617 = sin(qJ(1,3));
t622 = cos(qJ(3,3));
t623 = cos(qJ(1,3));
t640 = ((-t588 * t587 + t594 * t591) * t623 + (t587 * t591 + t594 * t588) * t617 + (t588 * t617 + t591 * t623) * m(3) * (-rSges(3,1) * t616 - rSges(3,2) * t622)) * t597;
t614 = legFrame(2,3);
t604 = sin(t614);
t607 = cos(t614);
t589 = -t604 * g(1) + t607 * g(2);
t592 = t607 * g(1) + t604 * g(2);
t595 = -qJ(2,2) * m(3) + (-rSges(2,3) - qJ(2,2)) * m(2) + t629;
t618 = sin(qJ(3,2));
t601 = t618 * pkin(3) + qJ(2,2);
t598 = 0.1e1 / t601;
t619 = sin(qJ(1,2));
t624 = cos(qJ(3,2));
t625 = cos(qJ(1,2));
t639 = ((-t589 * t587 + t595 * t592) * t625 + (t587 * t592 + t595 * t589) * t619 + (t589 * t619 + t592 * t625) * m(3) * (-rSges(3,1) * t618 - rSges(3,2) * t624)) * t598;
t615 = legFrame(1,3);
t605 = sin(t615);
t608 = cos(t615);
t590 = -t605 * g(1) + t608 * g(2);
t593 = t608 * g(1) + t605 * g(2);
t596 = -qJ(2,1) * m(3) + (-rSges(2,3) - qJ(2,1)) * m(2) + t629;
t620 = sin(qJ(3,1));
t602 = t620 * pkin(3) + qJ(2,1);
t599 = 0.1e1 / t602;
t621 = sin(qJ(1,1));
t626 = cos(qJ(3,1));
t627 = cos(qJ(1,1));
t638 = ((-t590 * t587 + t596 * t593) * t627 + (t587 * t593 + t596 * t590) * t621 + (t590 * t621 + t593 * t627) * m(3) * (-rSges(3,1) * t620 - rSges(3,2) * t626)) * t599;
t578 = -t588 * t623 + t591 * t617;
t637 = t578 * t597;
t579 = -t589 * t625 + t592 * t619;
t636 = t579 * t598;
t580 = -t590 * t627 + t593 * t621;
t635 = t580 * t599;
t630 = -m(2) - m(3);
t628 = rSges(3,1) * g(3);
t612 = 0.1e1 / t620;
t611 = 0.1e1 / t618;
t610 = 0.1e1 / t616;
t609 = pkin(6) + t642;
t586 = -t602 * t627 + t609 * t621;
t585 = -t601 * t625 + t609 * t619;
t584 = -t600 * t623 + t609 * t617;
t583 = t621 * t602 + t609 * t627;
t582 = t619 * t601 + t609 * t625;
t581 = t617 * t600 + t609 * t623;
t1 = [(-t605 * t621 + t608 * t627) * t638 + (-t604 * t619 + t607 * t625) * t639 + (-t603 * t617 + t606 * t623) * t640 - m(4) * g(1) + ((t583 * t608 - t605 * t586) * t635 + (t582 * t607 - t604 * t585) * t636 + (t581 * t606 - t603 * t584) * t637) * t630; (t605 * t627 + t608 * t621) * t638 + (t604 * t625 + t607 * t619) * t639 + (t603 * t623 + t606 * t617) * t640 - m(4) * g(2) + ((t605 * t583 + t586 * t608) * t635 + (t604 * t582 + t585 * t607) * t636 + (t603 * t581 + t584 * t606) * t637) * t630; -m(4) * g(3) + (t622 * t610 * t578 + t624 * t611 * t579 + t626 * t612 * t580) * t630 + (t612 * ((t580 * rSges(3,1) - t641) * t626 - t620 * (t580 * rSges(3,2) + t628)) + t611 * ((t579 * rSges(3,1) - t641) * t624 - t618 * (t579 * rSges(3,2) + t628)) + t610 * ((t578 * rSges(3,1) - t641) * t622 - t616 * (t578 * rSges(3,2) + t628))) / pkin(3) * m(3);];
taugX  = t1;
