% Calculate vector of inverse dynamics forces for parallel robot
% P6PRRRRR6V2G1P1A1
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% xDP [6x1]
%   Generalized platform velocities
% xDDP [6x1]
%   Generalized platform accelerations
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
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
% tauA [6x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in actuated joint coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-09-11 17:10
% Revision: 7993d029b5937b704dcf3fc7d8ae322038cdcbdd (2019-09-11)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauA = P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(10,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: xP has to be [6x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [6 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: xDP has to be [6x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [6 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: xDDP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: qJ has to be [3x6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PRRRRR6V2G1P1A1_invdyn_para_qa_slagn_vp1: Koppelpunkt has to be [6x3] (double)');

%% Function calls and calculation
tauX = P6PRRRRR6V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges);
Jinv = P6PRRRRR6V2G1P1A1_Jinv(xP, qJ, pkin, koppelP, legFrame);
tauA  = (Jinv') \ tauX;
tauA  = tauA;
